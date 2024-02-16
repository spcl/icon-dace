#include "cdi_fdb.h"

int cdi_fdb_dummy;

#ifdef HAVE_LIBFDB5

#include <string.h>
#include <stdlib.h>

#include "error.h"
#include "cdi_int.h"

void ensureBufferSize(size_t requiredSize, size_t *curSize, void **buffer);

void
decode_fdbitem(const char *fdbItem, KeyValueEntry *keyValue)
{
  keyValue->item = strdup(fdbItem);
  char *pItem = keyValue->item;
  int numKeys = 0;
  char **itemKeys = keyValue->keys;
  char **itemValues = keyValue->values;
  int len = strlen(pItem);
  int start = (*pItem == '{');
  itemKeys[0] = pItem + start;
  for (int i = start; i < len; i++)
    {
      if (pItem[i] == ',')
        {
          pItem[i] = 0;
          numKeys++;
          itemKeys[numKeys] = &pItem[++i];
        }
      else if (pItem[i] == '}')
        {
          pItem[i] = 0;
          numKeys++;
          if (pItem[i + 1] == '{')
            {
              itemKeys[numKeys] = &pItem[i + 2];
              i += 2;
            }
        }
      else if (i > start && i == (len - 1))
        {
          numKeys++;
        }
    }

  keyValue->numKeys = numKeys;
  // for (int i = 0; i < numKeys; i++)  printf("%d <%s>\n", i, itemKeys[i]);

  for (int i = 0; i < numKeys; i++)
    {
      char *itemKey = itemKeys[i];
      len = strlen(itemKey);
      for (int k = 0; k < len; k++)
        if (itemKey[k] == '=')
          {
            itemKey[k] = 0;
            itemValues[i] = &itemKey[k + 1];
            break;
          }

      // printf("key <%s> value <%s>\n", itemKeys[i], itemValues[i]);
    }
}

static int
fdb_request_add1(fdb_request_t *req, const char *param, const char *value)
{
  return fdb_request_add(req, param, &value, 1);
}

static void
fdbitem_to_request(const char *fdbItem, fdb_request_t *request)
{
  KeyValueEntry keyValue;
  keyValue.item = NULL;
  decode_fdbitem(fdbItem, &keyValue);

  for (int i = 0; i < keyValue.numKeys; i++)
    {
      // printf("key <%s> value <%s>\n", keyValue.keys[i], keyValue.values[i]);
      fdb_request_add1(request, keyValue.keys[i], keyValue.values[i]);
    }

  if (keyValue.item) free(keyValue.item);
}

int
fdb_fill_itemlist(fdb_handle_t *fdb, fdb_request_t *request, char ***itemList)
{
  const char **item = (const char **) malloc(sizeof(const char *));

  fdb_listiterator_t *it;
  fdb_new_listiterator(&it);

  fdb_list(fdb, request, it);

  int numItems = 0;
  while (true)
    {
      bool exist;
      fdb_listiterator_next(it, &exist, item);
      if (!exist) break;

      numItems++;
    }
  if (CDI_Debug) Message("numItems = %d", numItems);

  if (*itemList == NULL) *itemList = (char **) malloc(numItems * sizeof(char *));

  fdb_list(fdb, request, it);

  int itemNum = 0;
  while (true)
    {
      bool exist;
      fdb_listiterator_next(it, &exist, item);
      if (!exist) break;

      (*itemList)[itemNum++] = strdup(*item);
    }

  fdb_delete_listiterator(it);

  free(item);

  return numItems;
}

long
fdb_read_record(fdb_handle_t *fdb, char *item, size_t *buffersize, void **gribbuffer)
{
  // Message("%s", item);

  fdb_datareader_t *dataReader = NULL;
  fdb_new_datareader(&dataReader);
  fdb_request_t *singleRequest = NULL;
  fdb_new_request(&singleRequest);
  fdbitem_to_request(item, singleRequest);
  int status = fdb_retrieve(fdb, singleRequest, dataReader);
  fdb_delete_request(singleRequest);
  if (status != FDB_SUCCESS) Error("fdb_retrieve failed!");

  long recordSize = 0;
  fdb_datareader_open(dataReader, &recordSize);
  if (recordSize == 0) Error("fdb_datareader empty!");

  ensureBufferSize(recordSize, buffersize, gribbuffer);

  long readSize = 0;
  fdb_datareader_read(dataReader, *gribbuffer, recordSize, &readSize);
  // printf("fdb_datareader_read: size=%ld/%ld\n", recordSize, readSize);
  if (readSize != recordSize) Error("fdb_datareader_read failed!");

  fdb_datareader_close(dataReader);
  fdb_delete_datareader(dataReader);

  return recordSize;
}

static int
check_numKey(const char *key, int numKeys, int numItems)
{
  if (numKeys == 0)
    {
      Warning("Key %s is missing in all of the FDB records!", key);
      return -1;
    }
  else if (numKeys < numItems)
    {
      Warning("Key %s is missing in some of the FDB records!", key);
      return -2;
    }

  return 0;
}

int
check_keyvalueList(int numItems, KeyValueEntry *keyValueList)
{
  const char *searchKeys[] = { "date", "time", "param", "levtype" };
  const int numSearchKeys = sizeof(searchKeys) / sizeof(searchKeys[0]);
  int searchKeysCount[numSearchKeys];
  for (int k = 0; k < numSearchKeys; k++) searchKeysCount[k] = 0;

  for (int i = 0; i < numItems; i++)
    {
      int numKeys = keyValueList[i].numKeys;
      char **itemKeys = keyValueList[i].keys;
      for (int k = 0; k < numSearchKeys; k++)
        {

          for (int j = 0; j < numKeys; j++)
            {
              if (str_is_equal(itemKeys[j], searchKeys[k]))
                {
                  searchKeysCount[k]++;
                  break;
                }
            }
        }
    }

  int status = 0;
  for (int k = 0; k < numSearchKeys; k++)
    if (check_numKey(searchKeys[k], searchKeysCount[k], numItems) != 0) status = -1;

  return status;
}

typedef struct
{
  int date, time, param, levtype;
  int ilevel;
} CmpKeys;

void
record_info_entry_init(RecordInfoEntry *recordInfo)
{
  recordInfo->date = 0;
  recordInfo->time = 0;
  recordInfo->param = 0;
  recordInfo->levtype = 0;
  recordInfo->ilevel = 0;
}

static CmpKeys
set_cmpkeys(RecordInfoEntry *recordInfo)
{
  CmpKeys cmpKeys;
  cmpKeys.date = recordInfo->date;
  cmpKeys.time = recordInfo->time;
  cmpKeys.param = recordInfo->param;
  cmpKeys.levtype = recordInfo->levtype;
  cmpKeys.ilevel = recordInfo->ilevel;
  return cmpKeys;
}

static int
compare_cmpkeys(const CmpKeys *cmpKeys1, const CmpKeys *cmpKeys2)
{
  // clang-format off
  if (cmpKeys1->date == cmpKeys2->date &&
      cmpKeys1->time == cmpKeys2->time &&
      cmpKeys1->param == cmpKeys2->param &&
      cmpKeys1->levtype == cmpKeys2->levtype &&
      cmpKeys1->ilevel == cmpKeys2->ilevel)
    return 0;
  // clang-format on

  return -1;
}

int
get_num_records(int numItems, RecordInfoEntry *recordInfoList)
{
  const int date = recordInfoList[0].date;
  const int time = recordInfoList[0].time;

  int numRecords = 0;
  for (int i = 0; i < numItems; i++)
    {
      if (date == recordInfoList[i].date && time == recordInfoList[i].time)
        numRecords++;
      else
        break;
    }

  CmpKeys cmpKeys0 = set_cmpkeys(&recordInfoList[0]);
  for (int i = 1; i < numRecords; i++)
    {
      CmpKeys cmpKeys = set_cmpkeys(&recordInfoList[i]);
      if (compare_cmpkeys(&cmpKeys0, &cmpKeys) == 0)
        {
          numRecords = i;
          break;
        }
    }

  return numRecords;
}

enum
{
  levTypeUndef = 0,
  levTypeSFC,
  levTypeML,
  levTypePL
};

static int
get_ilevtype(const char *levtype)
{
  int ilevtype = levTypeUndef;

  // clang-format off
  if      (str_is_equal(levtype, "sfc")) ilevtype = levTypeSFC;
  else if (str_is_equal(levtype, "ml"))  ilevtype = levTypeML;
  else if (str_is_equal(levtype, "pl"))  ilevtype = levTypeML;
  // clang-format on

  return ilevtype;
}

void
decode_keyvalue(KeyValueEntry *keyValue, RecordInfoEntry *recordInfo)
{
  char **itemKeys = keyValue->keys;
  char **itemValues = keyValue->values;
  int numKeys = keyValue->numKeys;
  for (int i = 0; i < numKeys; i++)
    {
      // printf("key <%s> value <%s>\n", itemKeys[i], itemValues[i]);
      // clang-format off
      if      (str_is_equal(itemKeys[i], "date"))     recordInfo->date = atoi(itemValues[i]);
      else if (str_is_equal(itemKeys[i], "time"))     recordInfo->time = atoi(itemValues[i]);
      else if (str_is_equal(itemKeys[i], "param"))    recordInfo->param = atoi(itemValues[i]);
      else if (str_is_equal(itemKeys[i], "levtype"))  recordInfo->levtype = get_ilevtype(itemValues[i]);
      else if (str_is_equal(itemKeys[i], "levelist")) recordInfo->ilevel = atoi(itemValues[i]);
      // clang-format on
    }
}

int
remove_duplicate_timesteps(RecordInfoEntry *recordInfoList, int numRecords, int numTimesteps, int *timestepRecordOffset)
{
  int numTimestepsNew = numTimesteps;

  int date = recordInfoList[0].date;
  int time = recordInfoList[0].time;

  for (int i = 1; i < numTimesteps; ++i)
    {
      int k = 0;
      for (k = 0; k < numTimesteps; k++)
        {
          const int index = (i + k) * numRecords;
          if (date != recordInfoList[index].date || time != recordInfoList[index].time) break;
        }

      int index = i * numRecords;
      if (k > 0 && k < numTimesteps)
        {
          index = (i + k) * numRecords;
          int n = k;
          for (k = 0; k < n; k++)
            {
              Warning("Skip timestep %d", i + k + 1);
              numTimestepsNew--;
              for (int j = i; j < numTimestepsNew; j++) timestepRecordOffset[j] = timestepRecordOffset[j + 1];
            }
          i += k;
          if (i >= numTimesteps) break;
        }

      date = recordInfoList[index].date;
      time = recordInfoList[index].time;
    }

  return numTimestepsNew;
}

fdb_request_t *
create_fdb_request(const char *filename)
{
  size_t len = strlen(filename);
  if (len == 4) Error("Empty FDB request!");

  KeyValueEntry keyValue;
  keyValue.item = NULL;
  decode_fdbitem(filename + 4, &keyValue);

  if (keyValue.numKeys == 0) Error("Empty FDB request!");

  fdb_request_t *request = NULL;
  fdb_new_request(&request);

  bool classDefined = false;
  bool streamDefined = false;
  bool expverDefined = false;
  for (int i = 0; i < keyValue.numKeys; i++)
    {
      // clang-format off
      if      (str_is_equal(keyValue.keys[i], "class"))  classDefined = true;
      else if (str_is_equal(keyValue.keys[i], "stream")) streamDefined = true;
      else if (str_is_equal(keyValue.keys[i], "expver")) expverDefined = true;
      // clang-format on

      fdb_request_add1(request, keyValue.keys[i], keyValue.values[i]);
    }

  if (!classDefined) Error("FDB parameter <class> undefined!");
  if (!streamDefined) Error("FDB parameter <stream> undefined!");
  if (!expverDefined) Error("FDB parameter <expver> undefined!");

  /*
  fdb_request_add1(request, "class", "ea");
  fdb_request_add1(request, "expver", "0001");
  fdb_request_add1(request, "stream", "oper");
  fdb_request_add1(request, "domain", "g");
  fdb_request_add1(request, "date", "20180601");
  // fdb_request_add1(request, "time", "1800");
  fdb_request_add1(request, "type", "an");
  fdb_request_add1(request, "levtype", "sfc");
  fdb_request_add1(request, "step", "0");
  fdb_request_add1(request, "param", "139");
  // fdb_request_add1(request, "levelist", "300");
  */
  if (keyValue.item) free(keyValue.item);

  return request;
}

#endif
