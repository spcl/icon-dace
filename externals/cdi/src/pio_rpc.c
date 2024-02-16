#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdbool.h>

#include "cdi.h"
#include "dmemory.h"
#include "resource_handle.h"
/* FIXME: no longer needed when taxis updates are sent as separate data */
#include "taxis.h"

#include "pio_interface.h"
#include "pio_comm.h"
#include "pio_rpc.h"
#include "pio_util.h"

/* replace once size_t-based version of vlistInqVarSize is merged */
static size_t
cdiPioVlistInqVarSize(int vlistID, int varID)
{
  int zaxisID, gridID, tsteptype;
  vlistInqVar(vlistID, varID, &gridID, &zaxisID, &tsteptype);

  int nlevs = zaxisInqSize(zaxisID);

  int gridsize = gridInqSize(gridID);

  size_t size = (size_t) gridsize * (size_t) nlevs;

  return size;
}

struct clientBufSize
computeClientStreamBufSize(int streamID, const struct collSpec *collector)
{
  /* 1 record is filled in last to indicate number of records in total */
  struct clientBufSize rmaSizeSpec = { .bufSize = sizeof(struct winHeaderEntry), .numDataRecords = 1, .numRPCRecords = 0 };
  int vlistID = streamInqVlist(streamID);
  size_t nvars = (size_t) vlistNvars(vlistID);
  if (collector->partDesc)
    {
      /* the distribution of data is fully prescribed in this case,
       * i.e. memory needed can be computed exactly if the clients
       * also specified the element type as float or double, otherwise
       * plan for the double precision case */
      const Xt_idxlist *partDesc = collector->partDesc;
      const int *conversion = collector->conversion;
      for (size_t varID = 0; varID < nvars; ++varID)
        {
          size_t chunkSize = (size_t) (xt_idxlist_get_num_indices(partDesc[varID]));
          rmaSizeSpec.numDataRecords += 2;
          size_t elemSize = cdiPioElemSizeInference(varID, conversion);
          size_t chunkBytes = chunkSize * elemSize
                              /* re-align chunk to multiple of double size */
                              + elemSize
                              - 1
                              /* one header for data record, one for corresponding part
                               * descriptor*/
                              + 2 * sizeof(struct winHeaderEntry);
          rmaSizeSpec.bufSize += chunkBytes;
        }
    }
  else
    {
      /* the data distribution is unknown, use heuristic  */
      for (size_t varID = 0; varID < nvars; ++varID)
        {
          size_t chunkSize;
          {
            size_t varSize = cdiPioVlistInqVarSize(vlistID, (int) varID);
            chunkSize = (size_t) (ceilf(cdiPIOpartInflate_ * (float) varSize / (float) collector->numClients));
          }
          rmaSizeSpec.numDataRecords += 2;
          rmaSizeSpec.bufSize += chunkSize * sizeof(double)
                                 /* re-align chunk to multiple of double size */
                                 + sizeof(double)
                                 - 1
                                 /* one header for data record, one for corresponding part
                                  * descriptor*/
                                 + 2 * sizeof(struct winHeaderEntry)
                                 /* FIXME: heuristic for size of packed Xt_idxlist */
                                 + sizeof(Xt_int) * chunkSize * 3;
        }
    }

  // memory required for the function calls encoded
  // for remote execution
  // once per stream and timestep for each collector process
  // from one model process
  if (collector->sendRPCData)
    {
      int taxisID = vlistInqTaxis(vlistID);
      MPI_Comm comm = cdiPioInqInterComm();
      rmaSizeSpec.numRPCRecords = numRPCFuncs;
      rmaSizeSpec.bufSize += numRPCFuncs * sizeof(struct winHeaderEntry)
                             /* data part of streamDefTimestep */
                             + (size_t) (reshResourceGetPackSize(taxisID, &taxisOps, &comm));
    }
  rmaSizeSpec.bufSize = roundUpToMultiple(rmaSizeSpec.bufSize, PIO_WIN_ALIGN);
  return rmaSizeSpec;
}

void
cdiPioDestroyPartDescPreset(size_t nRanks, size_t nVars, struct partDescPreset *deco)
{
  Xt_idxlist(*restrict partDesc)[nRanks] = (Xt_idxlist(*)[nRanks]) deco->lists;
  Xt_uid(*restrict partDescUID)[nRanks] = (Xt_uid(*restrict)[nRanks]) deco->uids;
  deco->lists = NULL;
  deco->uids = NULL;
  for (size_t k = 0; k < nRanks; ++k)
    {
      for (size_t i = 0; i < nVars; ++i)
        {
          Xt_uid uid = partDescUID[i][k];
          if (uid != 0)
            {
              xt_idxlist_delete(partDesc[i][k]);
              for (size_t j = i + 1; j < nVars; ++j)
                if (partDescUID[j][k] == uid) partDescUID[j][k] = 0;
            }
        }
    }
  Free(deco->conversion);
  Free(partDesc);
  Free(partDescUID);
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
