#include "dmemory.h"
#include "stream_scan.h"

void
streamScanResizeRecords1(stream_t *streamptr)
{
  const int nrecords = streamptr->tsteps[0].nallrecs;
  if (nrecords < streamptr->tsteps[0].recordSize)
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records = (record_t *) Realloc(streamptr->tsteps[0].records, (size_t) nrecords * sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *) Malloc((size_t) nrecords * sizeof(int));
  streamptr->tsteps[0].nrecs = nrecords;
  for (int recID = 0; recID < nrecords; ++recID) streamptr->tsteps[0].recIDs[recID] = recID;
}

int
streamScanInitRecords2(stream_t *streamptr)
{
  const int nrecords = streamptr->tsteps[1].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) Malloc((size_t) nrecords * sizeof(int));
  streamptr->tsteps[1].nrecs = 0;

  for (int recID = 0; recID < nrecords; ++recID)
    {
      streamptr->tsteps[1].recIDs[recID] = -1;
      streamptr->tsteps[1].records[recID].position = streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[1].records[recID].size = streamptr->tsteps[0].records[recID].size;
    }

  return nrecords;
}

int
streamScanInitRecords(stream_t *streamptr, int tsID)
{
  const int nrecs = streamptr->tsteps[1].nrecs;

  streamptr->tsteps[tsID].nrecs = nrecs;
  streamptr->tsteps[tsID].recIDs = (int *) Malloc((size_t) nrecs * sizeof(int));

  for (int recID = 0; recID < nrecs; ++recID) streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

  return nrecs;
}

void
streamScanTimeConstAdjust(stream_t *streamptr, const taxis_t *taxis)
{
  const int vlistID = streamptr->vlistID;
  if (streamptr->ntsteps == 1 && cdiDateTime_isNull(taxis->vDateTime))
    {
      streamptr->ntsteps = 0;
      for (int varID = 0; varID < streamptr->nvars; ++varID) vlistDefVarTimetype(vlistID, varID, TIME_CONSTANT);
    }
}

void
streamScanTsFixNtsteps(stream_t *streamptr, off_t recpos)
{
  if (streamptr->ntsteps == -1)
    {
      const int tsID = tstepsNewEntry(streamptr);
      if (tsID != streamptr->rtsteps) Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID - 1].next = true;
      streamptr->tsteps[tsID].position = recpos;
    }
}
