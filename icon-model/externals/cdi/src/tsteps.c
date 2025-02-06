#include <limits.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"

static void
tstepsInitEntry(tsteps_t *tstep)
{
  tstep->recIDs = NULL;
  tstep->records = NULL;
  tstep->recordSize = 0;
  tstep->nrecs = 0;
  tstep->curRecID = CDI_UNDEFID;
  tstep->ncStepIndex = 0;
  tstep->position = 0;
  tstep->nallrecs = 0;
  tstep->next = 0;

  ptaxisInit(&(tstep->taxis));
}

int
tstepsNewEntry(stream_t *streamptr)
{
  const int tsID = streamptr->tstepsNextID++;
  int tstepsTableSize = streamptr->tstepsTableSize;
  tsteps_t *tstepsTable = streamptr->tsteps;

  // If the table overflows, double its size.
  if (tsID == tstepsTableSize)
    {
      if (tstepsTableSize == 0) tstepsTableSize = 1;
      if (tstepsTableSize <= INT_MAX / 2)
        tstepsTableSize *= 2;
      else if (tstepsTableSize < INT_MAX)
        tstepsTableSize = INT_MAX;
      else
        Error("Resizing of tstep table failed!");

      tstepsTable = (tsteps_t *) Realloc(tstepsTable, (size_t) tstepsTableSize * sizeof(tsteps_t));
    }

  streamptr->tstepsTableSize = tstepsTableSize;
  streamptr->tsteps = tstepsTable;

  tsteps_t *curTstep = &streamptr->tsteps[tsID];
  tstepsInitEntry(curTstep);

  return tsID;
}

void
cdi_create_timesteps(int numTimesteps, stream_t *streamptr)
{
  streamptr->ntsteps = (long) numTimesteps;
  if (numTimesteps < 0 || streamptr->tstepsTableSize > 0) return;

  int ntsteps = (numTimesteps == 0) ? 1 : numTimesteps;

  streamptr->tsteps = (tsteps_t *) Malloc((size_t) ntsteps * sizeof(tsteps_t));

  streamptr->tstepsTableSize = ntsteps;
  streamptr->tstepsNextID = ntsteps;

  for (int tsID = 0; tsID < ntsteps; tsID++)
    {
      tstepsInitEntry(&streamptr->tsteps[tsID]);
    }
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
