#ifndef CDI_PIO_IMPL_H
#define CDI_PIO_IMPL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "pio_conf.h"

enum IO_Server_command
{
  IO_Open_file,
  IO_Close_file,
  IO_Set_fp,
  IO_Send_buffer,
  IO_Finalize,
  IO_Sync_file,
  tagKey = 8, /* should be power of 2, must be
               * larger than IO_Finalize */
};

extern const char *const cdiPioCmdStrTab[];

struct syncMsg
{
  off_t amount;
  int fileID, command;
};

/* datatype to communicate an off_t */
extern MPI_Datatype cdiPioOffsetDt;
/* datatype to communicate a sync message */
extern MPI_Datatype cdiPioSyncMsgDt;

void cdiPioLookupOffsetDt(void);
void cdiPioCreateSyncMsgDt(void);
void cdiPioDestroySyncMsgDt(void);

struct fileOpTag
{
  int id;
  int command;
};

/* pio.c */
static inline int
encodeFileOpTag(int fileID, int command)
{
  return fileID * tagKey + command;
}

static inline struct fileOpTag
decodeFileOpTag(int tag)
{
  struct fileOpTag rtag = { .id = tag / tagKey, .command = tag % tagKey };
  return rtag;
}

/* pio_mpinonb.c */
void initMPINONB(void);

/* pio_mpi_fw_ordered.c */
void cdiPioFileWriteOrderedInit(void);

/* pio_mpi_fw_at_all.c */
void cdiPioFileWriteAtAllInit(void);

/* pio_mpi_fw_at_reblock.c */
void cdiPioFileWriteAtReblockInit(void);

/* common functionality for file split between collectors and writer(s) */
void pioSendInitialize(void);

/* pio_posixasynch.c */
#ifdef _POSIX_ASYNCHRONOUS_IO
void pioWriterAIO(void);
#endif

/* pio_posixfpguardsendrecv.c */
void initPOSIXFPGUARDSENDRECV(void);

/* pio_posixnonb.c */
void pioWriterStdIO(void);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
