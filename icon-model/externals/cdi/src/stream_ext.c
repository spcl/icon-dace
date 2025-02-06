#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_scan.h"
#include "stream_ext.h"
#include "varscan.h"
#include "extra.h"
#include "get_num_missvals.h"
#include "exse.h"

#ifdef HAVE_LIBEXTRA

static int
extInqDatatype(int prec, int number)
{
  if (number == 2)
    return (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
  else
    return (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
}

static void
extDefDatatype(int datatype, int *prec, int *number)
{
  if (datatype != CDI_DATATYPE_FLT32 && datatype != CDI_DATATYPE_FLT64 && datatype != CDI_DATATYPE_CPX32
      && datatype != CDI_DATATYPE_CPX64)
    datatype = CDI_DATATYPE_FLT32;

  *number = (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) ? 2 : 1;

  *prec = (datatype == CDI_DATATYPE_FLT64 || datatype == CDI_DATATYPE_CPX64) ? EXSE_DOUBLE_PRECISION : EXSE_SINGLE_PRECISION;
}

/* not used
int extInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int status;
  int fileID;
  int icode, ilevel;
  int zaxisID = -1;
  int header[4];
  int vlistID;
  void *extp = streamptr->record->objectp;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  *varID   = -1;
  *levelID = -1;

  status = extRead(fileID, extp);
  if ( status != 0 ) return 0;

  extInqHeader(extp, header);

  icode  = header[1];
  ilevel = header[2];

  *varID = vlistInqVarID(vlistID, icode);

  if ( *varID == CDI_UNDEFID ) Error("Code %d undefined", icode);

  zaxisID = vlistInqVarZaxis(vlistID, *varID);

  *levelID = zaxisInqLevelID(zaxisID, (double) ilevel);

  return 1;
}
*/

static void
ext_read_recordSP(stream_t *streamptr, float *data, size_t *numMissVals)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;

  int vrecID = streamptr->tsteps[tsID].curRecID;
  int recID = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *extp = streamptr->record->objectp;
  if (extRead(fileID, extp) < 0) Error("Failed to read EXTRA record");

  int header[4];
  extInqHeader(extp, header);
  extInqDataSP(extp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t size = gridInqSize(vlistInqVarGrid(vlistID, varID));

  if (vlistInqVarNumber(vlistID, varID) == CDI_REAL)
    *numMissVals = get_num_missvalsSP(size, data, missval);
  else
    *numMissVals = get_cplx_num_missvalsSP(size, data, missval);

  streamptr->numvals += size;
}

static void
ext_read_recordDP(stream_t *streamptr, double *data, size_t *numMissVals)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;

  int vrecID = streamptr->tsteps[tsID].curRecID;
  int recID = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *extp = streamptr->record->objectp;
  if (extRead(fileID, extp) < 0) Error("Failed to read EXTRA record");

  int header[4];
  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t size = gridInqSize(vlistInqVarGrid(vlistID, varID));

  if (vlistInqVarNumber(vlistID, varID) == CDI_REAL)
    *numMissVals = get_num_missvalsDP(size, data, missval);
  else
    *numMissVals = get_cplx_num_missvalsDP(size, data, missval);

  streamptr->numvals += size;
}

void
ext_read_record(stream_t *streamptr, int memtype, void *data, size_t *numMissVals)
{
  if (memtype == MEMTYPE_DOUBLE)
    ext_read_recordDP(streamptr, (double *) data, numMissVals);
  else
    ext_read_recordSP(streamptr, (float *) data, numMissVals);
}

void
extCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "EXTRA");
}

void
extDefRecord(stream_t *streamptr)
{
  Record *record = streamptr->record;

  int pdis, pcat, pnum;
  cdiDecodeParam(record->param, &pnum, &pcat, &pdis);

  int header[4];
  header[0] = record->vdate;
  header[1] = pnum;
  header[2] = record->ilevel;
  int gridID = record->gridID;
  cdi_check_gridsize_int_limit("EXTRA", gridInqSize(gridID));
  header[3] = (int) gridInqSize(gridID);

  extrec_t *extp = (extrec_t *) record->objectp;
  extDefDatatype(record->prec, &extp->prec, &extp->number);
  extDefHeader(extp, header);
}

static void
ext_write_recordSP(stream_t *streamptr, const float *data)
{
  void *extp = streamptr->record->objectp;
  extDefDataSP(extp, data);
  extWrite(streamptr->fileID, extp);
}

static void
ext_write_recordDP(stream_t *streamptr, const double *data)
{
  void *extp = streamptr->record->objectp;
  extDefDataDP(extp, data);
  extWrite(streamptr->fileID, extp);
}

void
ext_write_record(stream_t *streamptr, int memtype, const void *data)
{
  if (memtype == MEMTYPE_DOUBLE)
    ext_write_recordDP(streamptr, (const double *) data);
  else
    ext_write_recordSP(streamptr, (const float *) data);
}

static void
extAddRecord(stream_t *streamptr, int param, int level, size_t xysize, size_t recsize, off_t position, int prec, int number)
{
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int recID = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  record->size = recsize;
  record->position = position;
  record->param = param;
  record->ilevel = level;

  grid_t *grid = (grid_t *) Malloc(sizeof(*grid));
  grid_init(grid);
  cdiGridTypeInit(grid, GRID_GENERIC, xysize);
  grid->x.size = xysize;
  grid->y.size = 0;
  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 0);
  int gridID = gridAdded.Id;
  if (!gridAdded.isNew)
    {
      grid_free(grid);
      Free(grid);
    }

  int leveltype = ZAXIS_GENERIC;
  int datatype = extInqDatatype(prec, number);

  int varID, levelID = 0;
  varAddRecord(recID, param, gridID, leveltype, 0, level, 0, 0, 0, datatype, &varID, &levelID, TSTEP_INSTANT, 0, -1, NULL, NULL,
               NULL, NULL);

  record->varID = (short) varID;
  record->levelID = levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if (CDI_Debug) Message("varID = %d gridID = %d levelID = %d", varID, gridID, levelID);
}

static void
extScanTimestep1(stream_t *streamptr)
{
  CdiDateTime datetime0;
  cdiDateTime_init(&datetime0);
  int header[4];
  off_t recpos = 0;
  extrec_t *extp = (extrec_t *) streamptr->record->objectp;

  streamptr->curTsID = 0;

  int tsID = tstepsNewEntry(streamptr);
  if (tsID != 0) Error("Internal problem! tstepsNewEntry returns %d", tsID);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  int fileID = streamptr->fileID;

  int nrecs = 0;
  while (true)
    {
      recpos = fileGetPos(fileID);
      if (extRead(fileID, extp) != 0)
        {
          streamptr->ntsteps = 1;
          break;
        }

      size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

      extInqHeader(extp, header);

      int vdate = header[0];
      int vtime = 0;
      int rcode = header[1];
      int rlevel = header[2];
      int rxysize = header[3];
      int param = cdiEncodeParam(rcode, 255, 255);
      CdiDateTime datetime = cdiDateTime_set(vdate, vtime);

      if (nrecs == 0)
        {
          datetime0 = datetime;
          taxis->vDateTime = datetime;
        }
      else
        {
          record_t *records = streamptr->tsteps[tsID].records;
          for (int recID = 0; recID < nrecs; recID++)
            if (param == records[recID].param && rlevel == records[recID].ilevel) goto tstepScanLoopFinished;

          if (cdiDateTime_isNE(datetime, datetime0)) Warning("Inconsistent verification time for code %d level %d", rcode, rlevel);
        }

      nrecs++;

      if (CDI_Debug) Message("%4d%8d%4d%8d%8d%6d", nrecs, (int) recpos, rcode, rlevel, vdate, vtime);

      extAddRecord(streamptr, param, rlevel, rxysize, recsize, recpos, extp->prec, extp->number);
    }

tstepScanLoopFinished:
  streamptr->rtsteps = 1;

  cdi_generate_vars(streamptr);

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type = TAXIS_ABSOLUTE;
  taxis->rDateTime = taxis->vDateTime;

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  vlist_check_contents(vlistID);

  streamScanResizeRecords1(streamptr);

  streamScanTsFixNtsteps(streamptr, recpos);
  streamScanTimeConstAdjust(streamptr, taxis);
}

static int
extScanTimestep2(stream_t *streamptr)
{
  int header[4];
  off_t recpos = 0;
  void *extp = streamptr->record->objectp;

  streamptr->curTsID = 1;

  int fileID = streamptr->fileID;
  int vlistID = streamptr->vlistID;

  int tsID = streamptr->rtsteps;
  if (tsID != 1) Error("Internal problem! unexpected timestep %d", tsID + 1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);
  record_t *records = streamptr->tsteps[tsID].records;

  int nrecords = streamScanInitRecords2(streamptr);

  for (int rindex = 0; rindex <= nrecords; rindex++)
    {
      recpos = fileGetPos(fileID);
      if (extRead(fileID, extp) != 0)
        {
          streamptr->ntsteps = 2;
          break;
        }

      size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

      extInqHeader(extp, header);

      int vdate = header[0];
      int vtime = 0;
      int rcode = header[1];
      int rlevel = header[2];
      int param = cdiEncodeParam(rcode, 255, 255);

      if (rindex == 0)
        {
          taxis->type = TAXIS_ABSOLUTE;
          taxis->vDateTime = cdiDateTime_set(vdate, vtime);
        }

      bool nextstep = false;
      int recID;
      for (recID = 0; recID < nrecords; recID++)
        {
          if (param == records[recID].param && rlevel == records[recID].ilevel)
            {
              if (records[recID].used)
                {
                  nextstep = true;
                }
              else
                {
                  records[recID].used = true;
                  streamptr->tsteps[tsID].recIDs[rindex] = recID;
                }
              break;
            }
        }
      if (recID == nrecords)
        {
          Warning("Code %d level %d not found at timestep %d", rcode, rlevel, tsID + 1);
          return CDI_EUFSTRUCT;
        }

      if (nextstep) break;

      if (CDI_Debug) Message("%4d%8d%4d%8d%8d%6d", rindex + 1, (int) recpos, rcode, rlevel, vdate, vtime);

      if (param != records[recID].param || rlevel != records[recID].ilevel)
        {
          Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                  records[recID].ilevel, rlevel);
          return CDI_EUFSTRUCT;
        }

      records[recID].position = recpos;
      records[recID].size = recsize;
    }

  int nrecs = 0;
  for (int recID = 0; recID < nrecords; recID++)
    {
      if (records[recID].used)
        nrecs++;
      else
        vlistDefVarTimetype(vlistID, records[recID].varID, TIME_CONSTANT);
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  streamScanTsFixNtsteps(streamptr, recpos);

  return 0;
}

int
extInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;

  extScanTimestep1(streamptr);

  int status = (streamptr->ntsteps == -1) ? extScanTimestep2(streamptr) : 0;

  fileSetPos(streamptr->fileID, 0, SEEK_SET);

  return status;
}

static long
extScanTimestep(stream_t *streamptr)
{
  int header[4];
  off_t recpos = 0;
  int nrecs = 0;
  void *extp = streamptr->record->objectp;

  int tsID = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if (streamptr->tsteps[tsID].recordSize == 0)
    {
      cdi_create_records(streamptr, tsID);
      record_t *records = streamptr->tsteps[tsID].records;

      nrecs = streamScanInitRecords(streamptr, tsID);

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for (int rindex = 0; rindex <= nrecs; rindex++)
        {
          recpos = fileGetPos(fileID);
          if (extRead(fileID, extp) != 0)
            {
              streamptr->ntsteps = streamptr->rtsteps + 1;
              break;
            }

          size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

          extInqHeader(extp, header);

          int vdate = header[0];
          int vtime = 0;
          int rcode = header[1];
          int rlevel = header[2];
          int param = cdiEncodeParam(rcode, 255, 255);

          // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
          if (rindex == nrecs) continue;
          int recID = streamptr->tsteps[tsID].recIDs[rindex];

          if (rindex == 0)
            {
              taxis->type = TAXIS_ABSOLUTE;
              taxis->vDateTime = cdiDateTime_set(vdate, vtime);
            }

          if (param != records[recID].param || rlevel != records[recID].ilevel)
            {
              Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                      records[recID].ilevel, rlevel);
              Error("Invalid, unsupported or inconsistent record structure!");
            }

          records[recID].position = recpos;
          records[recID].size = recsize;

          if (CDI_Debug) Message("%4d%8d%4d%8d%8d%6d", rindex, (int) recpos, rcode, rlevel, vdate, vtime);
        }

      streamptr->rtsteps++;

      if (streamptr->ntsteps != streamptr->rtsteps)
        {
          tsID = tstepsNewEntry(streamptr);
          if (tsID != streamptr->rtsteps) Error("Internal error. tsID = %d", tsID);

          streamptr->tsteps[tsID - 1].next = true;
          streamptr->tsteps[tsID].position = recpos;
        }

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;
    }

  if (nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs)
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return streamptr->ntsteps;
}

int
extInqTimestep(stream_t *streamptr, int tsID)
{
  if (tsID == 0 && streamptr->rtsteps == 0) Error("Call to cdiInqContents missing!");

  if (CDI_Debug) Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  long ntsteps = CDI_UNDEFID;
  while ((tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID) ntsteps = extScanTimestep(streamptr);

  int nrecs = 0;
  if (!(tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID))
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}

void
extReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, size_t *numMissVals)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  void *extp = streamptr->record->objectp;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  int tsid = streamptr->curTsID;

  off_t currentfilepos = fileGetPos(fileID);

  int recID = streamptr->vars[varID].recordTable[0].recordID[levID];
  off_t recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  extRead(fileID, extp);
  int header[4];
  extInqHeader(extp, header);
  extInqDataDP(extp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  if (vlistInqVarNumber(vlistID, varID) == CDI_REAL)
    *numMissVals = get_num_missvalsDP(gridsize, data, missval);
  else
    *numMissVals = get_cplx_num_missvalsDP(gridsize, data, missval);
}

void
extReadVarDP(stream_t *streamptr, int varID, double *data, size_t *numMissVals)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs = (size_t) streamptr->vars[varID].recordTable[0].nlevs;

  for (size_t levID = 0; levID < nlevs; levID++)
    extReadVarSliceDP(streamptr, varID, (int) levID, &data[levID * gridsize], numMissVals);
}

void
extWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;
  CdiDateTime vDateTime = streamptr->tsteps[tsID].taxis.vDateTime;

  int pdis, pcat, pnum;
  cdiDecodeParam(vlistInqVarParam(vlistID, varID), &pnum, &pcat, &pdis);

  int header[4];
  header[0] = (int) cdiDate_get(vDateTime.date);
  header[1] = pnum;
  header[2] = (int) lround(zaxisInqLevel(vlistInqVarZaxis(vlistID, varID), levID));
  int gridID = vlistInqVarGrid(vlistID, varID);
  cdi_check_gridsize_int_limit("EXTRA", gridInqSize(gridID));
  header[3] = (int) gridInqSize(gridID);

  extrec_t *extp = (extrec_t *) streamptr->record->objectp;
  extDefDatatype(vlistInqVarDatatype(vlistID, varID), &extp->prec, &extp->number);
  extDefHeader(extp, header);

  extDefDataDP(extp, data);
  extWrite(fileID, extp);
}

void
extWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs = (size_t) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));

  for (size_t levID = 0; levID < nlevs; levID++) extWriteVarSliceDP(streamptr, varID, (int) levID, &data[levID * gridsize]);
}

#endif /* HAVE_LIBEXTRA */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
