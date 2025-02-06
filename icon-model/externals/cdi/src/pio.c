#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "cdi_int.h"
#include "pio.h"
#include "cdipio.h"
#include "pio_impl.h"
#include "pio_util.h"

const char *const cdiPioCmdStrTab[] = {
  "IO_Open_file", "IO_Close_file", "IO_Set_fp", "IO_Send_buffer", "IO_Finalize", "IO_Sync_file",
};

int cdiPioExtraNSKeys[cdiPioExtraNSKeysSize];

MPI_Datatype cdiPioOffsetDt = MPI_DATATYPE_NULL;
MPI_Datatype cdiPioSyncMsgDt = MPI_DATATYPE_NULL;

void
cdiPioLookupOffsetDt(void)
{
  xmpi(MPI_Type_match_size(MPI_TYPECLASS_INTEGER, (int) (sizeof(off_t)), &cdiPioOffsetDt));
}

void
cdiPioCreateSyncMsgDt(void)
{
  if (cdiPioOffsetDt == MPI_DATATYPE_NULL) cdiPioLookupOffsetDt();
  struct syncMsg dummy;
  int bl[2] = { 1, 2 };
  MPI_Aint displ[2];
  displ[0] = 0;
  displ[1] = (unsigned char *) &dummy.fileID - (unsigned char *) &dummy.amount;
  MPI_Datatype dt[2] = { cdiPioOffsetDt, MPI_INT };
  xmpi(MPI_Type_create_struct(2, bl, displ, dt, &cdiPioSyncMsgDt));
  xmpi(MPI_Type_commit(&cdiPioSyncMsgDt));
}

void
cdiPioDestroySyncMsgDt(void)
{
  MPI_Type_free(&cdiPioSyncMsgDt);
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
