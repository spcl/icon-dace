#include "dmemory.h"
#include "cdi_int.h"

static void
cdf_init_timestep(tsteps_t *timeStep, int numRecs, int numRecsAvail)
{
  timeStep->records = (record_t *) Malloc((size_t) numRecs * sizeof(record_t));
  timeStep->nrecs = numRecsAvail;
  timeStep->nallrecs = numRecs;
  timeStep->recordSize = numRecs;
  timeStep->curRecID = CDI_UNDEFID;
}

static int
cdf_get_numRecsAvail(int vlistID)
{
  int numRecsAvail = 0;
  int numVars = vlistNvars(vlistID);
  for (int varID = 0; varID < numVars; varID++)
    {
      if (vlistInqVarTimetype(vlistID, varID) != TIME_CONSTANT)
        {
          numRecsAvail += zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
        }
    }
  return numRecsAvail;
}

static void
cdf_init_records_step0(int numRecs, int *recIDs, record_t *records, int vlistID)
{
  for (int recID = 0; recID < numRecs; recID++) recIDs[recID] = recID;

  int numVars = vlistNvars(vlistID);
  for (int varID = 0, recID = 0; varID < numVars; varID++)
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int nlevels = zaxisInqSize(zaxisID);
      for (int levelID = 0; levelID < nlevels; levelID++)
        {
          recordInitEntry(&records[recID]);
          records[recID].varID = (short) varID;
          records[recID].levelID = levelID;
          recID++;
        }
    }
}

static void
cdf_init_records_step1(int numRecs, int *recIDs, record_t *records, int vlistID)
{
  for (int recID = 0, vrecID = 0; recID < numRecs; recID++)
    {
      if (vlistInqVarTimetype(vlistID, records[recID].varID) != TIME_CONSTANT)
        {
          recIDs[vrecID++] = recID;
        }
    }
}

void
cdf_create_records(stream_t *streamptr, int tsID)
{
  if (tsID < 0 || (tsID >= streamptr->ntsteps && tsID > 0)) return;

  if (streamptr->tsteps[tsID].nallrecs > 0) return;

  int vlistID = streamptr->vlistID;

  tsteps_t *sourceTstep = streamptr->tsteps;
  tsteps_t *destTstep = sourceTstep + tsID;

  int numRecs = vlistNrecs(vlistID);
  if (numRecs <= 0) return;

  if (tsID == 0)
    {
      int numRecsAvail = numRecs;  // use all records at first timestep

      streamptr->nrecs += numRecs;

      cdf_init_timestep(destTstep, numRecs, numRecsAvail);

      destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
      cdf_init_records_step0(numRecs, destTstep->recIDs, destTstep->records, vlistID);
    }
  else if (tsID == 1)
    {
      int numRecsAvail = cdf_get_numRecsAvail(vlistID);

      streamptr->nrecs += numRecsAvail;

      cdf_init_timestep(destTstep, numRecs, numRecsAvail);

      memcpy(destTstep->records, sourceTstep->records, (size_t) numRecs * sizeof(record_t));

      if (numRecsAvail)
        {
          destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
          cdf_init_records_step1(numRecs, destTstep->recIDs, destTstep->records, vlistID);
        }
    }
  else
    {
      if (streamptr->tsteps[1].records == 0) cdf_create_records(streamptr, 1);

      int numRecsAvail = streamptr->tsteps[1].nrecs;

      streamptr->nrecs += numRecsAvail;

      cdf_init_timestep(destTstep, numRecs, numRecsAvail);

      memcpy(destTstep->records, sourceTstep->records, (size_t) numRecs * sizeof(record_t));

      if (numRecsAvail)
        {
          destTstep->recIDs = (int *) Malloc((size_t) numRecsAvail * sizeof(int));
          memcpy(destTstep->recIDs, streamptr->tsteps[1].recIDs, (size_t) numRecsAvail * sizeof(int));
        }
    }
}
