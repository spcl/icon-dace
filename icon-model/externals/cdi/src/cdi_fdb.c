#include "cdi_fdb.h"

int cdi_fdb_dummy;

#ifdef HAVE_LIBFDB5

#include <string.h>
#include <stdlib.h>

#include "error.h"
#include "cdi_int.h"

void ensureBufferSize(size_t requiredSize, size_t *curSize, void **buffer);

void
check_fdb_error(int errorNum)
{
  // This routine provides a simple interface to FDB error message routine.
  if (errorNum != FDB_SUCCESS) Error("%s", fdb_error_string(errorNum));
}

void
cdi_fdb_delete_kvlist(int numItems, fdbKeyValueEntry *keyValueList)
{
  for (int i = 0; i < numItems; ++i)
    {
      for (int k = 0; k < keyValueList[i].numKeys; ++k)
        {
          if (keyValueList[i].keys[k]) free(keyValueList[i].keys[k]);
          if (keyValueList[i].values[k]) free(keyValueList[i].values[k]);
        }
    }
  free(keyValueList);
}

void
decode_fdbitem(const char *fdbItem, KeyValueItem *keyValue)
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

  // keyValue->numKeys = numKeys;
  // for (int i = 0; i < numKeys; i++)  printf("%d <%s>\n", i, itemKeys[i]);

  int n = 0;
  for (int i = 0; i < numKeys; i++)
    {
      char *itemKey = itemKeys[i];
      len = strlen(itemKey);
      for (int k = 0; k < len; k++)
        if (itemKey[k] == '=')
          {
            n++;
            itemKey[k] = 0;
            itemValues[i] = &itemKey[k + 1];
            break;
          }

      // printf("key <%s> value <%s>\n", itemKeys[i], itemValues[i]);
    }

  keyValue->numKeys = n;
}

static int
fdb_request_add1(fdb_request_t *req, const char *param, const char *value)
{
  return fdb_request_add(req, param, &value, 1);
}

static void
fdbitem_to_request(const fdbKeyValueEntry *keyValue, fdb_request_t *request)
{
  // printf("numKeys: %d\n", keyValue->numKeys);
  for (int i = 0; i < keyValue->numKeys; i++)
    {
      // printf("%d: key <%s> value <%s>\n", i, keyValue->keys[i], keyValue->values[i]);
      check_fdb_error(fdb_request_add1(request, keyValue->keys[i], keyValue->values[i]));
    }
}

int
cdi_fdb_fill_kvlist(fdb_handle_t *fdb, fdb_request_t *request, fdbKeyValueEntry **pkeyValueList)
{
  fdb_listiterator_t *it;

  bool reportDuplicatedListElements = false;
  check_fdb_error(fdb_list(fdb, request, &it, reportDuplicatedListElements));

  int numItems = 0;
  while (true)
    {
      int err = fdb_listiterator_next(it);
      if (err != FDB_SUCCESS) break;

      numItems++;
    }
  if (CDI_Debug) Message("numItems = %d", numItems);
  // Message("numItems = %d", numItems);

  if (*pkeyValueList == NULL) *pkeyValueList = (fdbKeyValueEntry *) malloc(numItems * sizeof(fdbKeyValueEntry));

  check_fdb_error(fdb_list(fdb, request, &it, reportDuplicatedListElements));

  int itemNum = 0;
  while (true)
    {
      int err = fdb_listiterator_next(it);
      if (err != FDB_SUCCESS) break;

      // const char *uri;
      // size_t off, attr_len;
      // fdb_listiterator_attrs(it, &uri, &off, &attr_len);
      // printf("uri=%s, off=%zu, attr_len=%zu\n", uri, off, attr_len);

      fdb_split_key_t *sk = NULL;
      check_fdb_error(fdb_new_splitkey(&sk));
      check_fdb_error(fdb_listiterator_splitkey(it, sk));

      int keyNum = 0;
      while (true)
        {
          const char *k;
          const char *v;
          size_t l;
          bool checkLevel = true;
          err = fdb_splitkey_next_metadata(sk, &k, &v, checkLevel ? &l : NULL);
          if (err != FDB_SUCCESS) break;
          // printf("k, v, l %s %s %zu\n", k, v, l);
          if (keyNum < 32)
            {
              (*pkeyValueList)[itemNum].keys[keyNum] = strdup(k);
              (*pkeyValueList)[itemNum].values[keyNum] = strdup(v);
              keyNum++;
            }
        }

      (*pkeyValueList)[itemNum++].numKeys = keyNum;

      check_fdb_error(fdb_delete_splitkey(sk));
    }

  check_fdb_error(fdb_delete_listiterator(it));

  return numItems;
}

long
cdi_fdb_read_record(fdb_handle_t *fdb, const fdbKeyValueEntry *keyValue, size_t *buffersize, void **gribbuffer)
{
  fdb_datareader_t *dataReader = NULL;
  check_fdb_error(fdb_new_datareader(&dataReader));
  fdb_request_t *singleRequest = NULL;
  check_fdb_error(fdb_new_request(&singleRequest));
  fdbitem_to_request(keyValue, singleRequest);
  check_fdb_error(fdb_retrieve(fdb, singleRequest, dataReader));
  check_fdb_error(fdb_delete_request(singleRequest));

  long recordSize = 0;
  check_fdb_error(fdb_datareader_open(dataReader, &recordSize));
  if (recordSize == 0) Error("fdb_datareader empty!");

  ensureBufferSize(recordSize, buffersize, gribbuffer);

  long readSize = 0;
  check_fdb_error(fdb_datareader_read(dataReader, *gribbuffer, recordSize, &readSize));
  // printf("fdb_datareader_read: size=%ld/%ld\n", recordSize, readSize);
  if (readSize != recordSize) Error("fdb_datareader_read failed!");

  check_fdb_error(fdb_datareader_close(dataReader));
  check_fdb_error(fdb_delete_datareader(dataReader));

  return recordSize;
}

static int
check_numKey(const char *key, int numKeys, int numItems)
{
  if (numKeys == 0)
    {
      Warning("Key %s is missing in all FDB records!", key);
      return -1;
    }
  else if (numKeys < numItems)
    {
      Warning("Key %s is missing in some FDB records!", key);
      return -2;
    }

  return 0;
}
/*
int
check_keyvalueList(int numItems, fdbKeyValueEntry *keyValueList)
{
  const char *searchKeys[] = { "param", "levtype", "date", "time" };
  int numSearchKeys = sizeof(searchKeys) / sizeof(searchKeys[0]);
  int searchKeysCount[numSearchKeys];
  for (int k = 0; k < numSearchKeys; k++) searchKeysCount[k] = 0;

  // for (int i = 0; i < numItems; i++)
  for (int i = 0; i < 1; i++)
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
*/
void
print_keyvalueList(int numItems, fdbKeyValueEntry *keyValueList)
{
  for (int i = 0; i < numItems; ++i)
    {
      printf("item=%d ", i + 1);
      const fdbKeyValueEntry *e = &keyValueList[i];
      for (int k = 0; k < e->numKeys; ++k) printf("%s%s=%s", (k > 0) ? "," : "", e->keys[k], e->values[k]);
      printf("\n");
    }
}

void
print_keyvalueList_sorted(int numItems, fdbKeyValueEntry *keyValueList, RecordInfoEntry *recordInfoList)
{
  for (int i = 0; i < numItems; ++i)
    {
      int fdbIndex = recordInfoList[i].fdbIndex;
      printf("item=%d ", fdbIndex + 1);
      const fdbKeyValueEntry *e = &keyValueList[fdbIndex];
      for (int k = 0; k < e->numKeys; ++k) printf("%s%s=%s", (k > 0) ? "," : "", e->keys[k], e->values[k]);
      printf("\n");
    }
}

static int
cmp_datetime(const void *e1, const void *e2)
{
  const RecordInfoEntry *x = (const RecordInfoEntry *) e1, *y = (const RecordInfoEntry *) e2;
  int64_t datetime1 = (int64_t) x->date * 100000 + x->time;
  int64_t datetime2 = (int64_t) y->date * 100000 + y->time;

  if (datetime1 < datetime2) return -1;
  if (datetime1 > datetime2) return 1;
  return 0;
}

static bool
isSorted_dateTime(int numItems, RecordInfoEntry *recordInfo)
{
  int64_t datetime1 = (int64_t) recordInfo[0].date * 100000 + recordInfo[0].time;
  for (int i = 1; i < numItems; ++i)
    {
      int64_t datetime2 = (int64_t) recordInfo[i].date * 100000 + recordInfo[i].time;
      if (datetime1 > datetime2) return false;
      datetime1 = datetime2;
    }

  return true;
}

void
cdi_fdb_sort_datetime(int numItems, RecordInfoEntry *recordInfo)
{
  if (!isSorted_dateTime(numItems, recordInfo)) qsort(recordInfo, numItems, sizeof(recordInfo[0]), cmp_datetime);
}

static void
record_info_entry_init(RecordInfoEntry *recordInfo)
{
  recordInfo->fdbIndex = -1;
  recordInfo->date = 0;
  recordInfo->time = 0;
  recordInfo->param = 0;
  recordInfo->levtype = 0;
  recordInfo->ilevel = 0;
}
/*
static int
compare_record_info_entry(const RecordInfoEntry *r1, const RecordInfoEntry *r2)
{
  // clang-format off
  if (r1->date    == r2->date &&
      r1->time    == r2->time &&
      r1->param   == r2->param &&
      r1->levtype == r2->levtype &&
      r1->ilevel  == r2->ilevel)
    return 0;
  // clang-format on

  return -1;
}
*/
int
get_num_records(int numItems, RecordInfoEntry *recordInfo)
{
  int numRecords = 0;
  for (int i = 0; i < numItems; i++)
    {
      if (recordInfo[0].date != recordInfo[i].date || recordInfo[0].time != recordInfo[i].time) break;
      numRecords++;
    }

  int numTimesteps = numItems / numRecords;
  if (numTimesteps * numRecords != numItems) return 0;

  for (int k = 1; k < numTimesteps; ++k)
    {
      int date = recordInfo[k * numRecords].date;
      int time = recordInfo[k * numRecords].time;
      for (int i = 1; i < numRecords; i++)
        {
          int index = k * numRecords + i;
          if (date != recordInfo[index].date || time != recordInfo[index].time) return 0;
        }
    }
  /*
  for (int i = 1; i < numRecords; i++)
    {
      if (compare_record_info_entry(&recordInfo[0], &recordInfo[i]) == 0)
        {
          numRecords = i;
          break;
        }
    }
  */
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
  else if (str_is_equal(levtype, "pl"))  ilevtype = levTypePL;
  // clang-format on

  return ilevtype;
}

int
decode_keyvalue(int numItems, fdbKeyValueEntry *keyValueList, RecordInfoEntry *recordInfoList)
{
  int numKeyDate = 0;
  int numKeyTime = 0;
  int numKeyParam = 0;
  int numKeyLtype = 0;
  for (int i = 0; i < numItems; ++i)
    {
      fdbKeyValueEntry *keyValue = &keyValueList[i];
      RecordInfoEntry *rentry = &recordInfoList[i];
      record_info_entry_init(rentry);
      rentry->fdbIndex = i;
      char **keys = keyValue->keys;
      char **values = keyValue->values;
      bool foundDate = false;
      bool foundTime = false;
      bool foundParam = false;
      bool foundLtype = false;
      bool foundLlist = false;
      for (int i = 0; i < keyValue->numKeys; i++)
        {
          // printf("key <%s> value <%s>\n", itemKeys[i], itemValues[i]);
          // clang-format off
          if      (!foundDate  && str_is_equal(keys[i], "date"))     { foundDate = true;  numKeyDate++;  rentry->date = atoi(values[i]); }
          else if (!foundTime  && str_is_equal(keys[i], "time"))     { foundTime = true;  numKeyTime++;  rentry->time = atoi(values[i]); }
          else if (!foundParam && str_is_equal(keys[i], "param"))    { foundParam = true; numKeyParam++; rentry->param = atoi(values[i]); }
          else if (!foundLtype && str_is_equal(keys[i], "levtype"))  { foundLtype = true; numKeyLtype++; rentry->levtype = get_ilevtype(values[i]); }
          else if (!foundLlist && str_is_equal(keys[i], "levelist")) { foundLlist = true; rentry->ilevel = atoi(values[i]); }
          // clang-format on
          if (foundDate && foundTime && foundParam && foundLtype && foundLlist) break;
        }
    }

  int status = 0;
  if (check_numKey("date", numKeyDate, numItems) != 0) status = -1;
  if (check_numKey("time", numKeyTime, numItems) != 0) status = -1;
  if (check_numKey("param", numKeyParam, numItems) != 0) status = -1;
  if (check_numKey("levtype", numKeyLtype, numItems) != 0) status = -1;

  return status;
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
          int index = (i + k) * numRecords;
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
cdi_create_fdb_request(const char *filename)
{
  size_t len = strlen(filename);
  if (len == 6) Error("Empty FDB request!");

  KeyValueItem keyValue;
  keyValue.item = NULL;
  decode_fdbitem(filename + 6, &keyValue);

  if (keyValue.numKeys == 0) Error("Empty FDB request!");

  fdb_request_t *request = NULL;
  fdb_new_request(&request);

  bool expverDefined = false;
  bool classDefined = false;
  for (int k = 0; k < keyValue.numKeys; k++)
    {
      // clang-format off
      if      (!expverDefined && str_is_equal(keyValue.keys[k], "expver")) expverDefined = true;
      else if (!classDefined  && str_is_equal(keyValue.keys[k], "class")) classDefined = true;
      // clang-format on

      check_fdb_error(fdb_request_add1(request, keyValue.keys[k], keyValue.values[k]));
    }

  if (!expverDefined) Error("FDB parameter <expver> undefined!");
  if (!classDefined) Error("FDB parameter <class> undefined!");

  if (keyValue.item) free(keyValue.item);

  return request;
}

#endif
