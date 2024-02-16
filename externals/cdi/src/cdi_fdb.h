#ifndef CDI_FDB_H
#define CDI_FDB_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern int cdi_fdb_dummy;

#ifdef HAVE_LIBFDB5

#include <fdb5/api/fdb_c.h>

typedef struct
{
  char *item;
  char *keys[32];
  char *values[32];
  int numKeys;
} KeyValueEntry;

typedef struct
{
  int date;
  int time;
  int param;
  int levtype;
  int ilevel;
} RecordInfoEntry;

void decode_fdbitem(const char *fdbItem, KeyValueEntry *keyValue);
int fdb_fill_itemlist(fdb_handle_t *fdb, fdb_request_t *request, char ***itemList);
long fdb_read_record(fdb_handle_t *fdb, char *item, size_t *buffersize, void **gribbuffer);
int check_keyvalueList(int numItems, KeyValueEntry *keyValueList);
void record_info_entry_init(RecordInfoEntry *recordInfo);
int get_num_records(int numItems, RecordInfoEntry *recordInfoList);
void decode_keyvalue(KeyValueEntry *keyValue, RecordInfoEntry *recordInfo);
int remove_duplicate_timesteps(RecordInfoEntry *recordInfoList, int numRecords, int numTimesteps, int *timestepRecordOffset);
fdb_request_t *create_fdb_request(const char *filename);

#endif

#endif /* CDI_FDB_H */
