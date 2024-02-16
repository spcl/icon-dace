#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "cdi.h"
#include "cdi_uuid.h"
#include "dmemory.h"
#include "error.h"
#include "resource_handle.h"
#include "resource_unpack.h"

#ifdef USE_MPI
#include <mpi.h>
#include "cdipio.h"
#include "pio_serialize.h"
#include "pio_util.h"
#else
typedef int MPI_Comm;
#endif

enum
{
  DOUBLE_PRECISION = CDI_DATATYPE_FLT64,
  nlon = 12,
  nlat = 6,
  nlev = 5,
  ntsteps = 3
};

static double lons[nlon] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
static double lats[nlat] = { -75, -45, -15, 15, 45, 75 };
static double levs[nlev] = { 101300, 92500, 85000, 50000, 20000 };

static int
defineGrid(void)
{
  int mask_vec[nlon * nlat];
  const int *mp = &mask_vec[0];
  double area_vec[nlon * nlat];
  const double *ap = &area_vec[0];
  int i;

  int gridID = gridCreate(GRID_LONLAT, nlon * nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);
  gridDefNvertex(gridID, 1);
  gridDefXbounds(gridID, lons);
  gridDefYbounds(gridID, lats);
  for (i = 0; i < nlon * nlat; i++) mask_vec[i] = i % 2;
  gridDefMaskGME(gridID, mp);
  for (i = 0; i < nlon * nlat; i++) mask_vec[i] = 1;
  gridDefMask(gridID, mp);

  cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, "myXname");
  cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, "myYname");
  cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, "myXlongname");
  cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, "myYlongname");
  cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, "myXunits");
  cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, "myYunits");

  gridDefDatatype(gridID, DOUBLE_PRECISION);
  gridDefTrunc(gridID, 1);
  gridDefParamGME(gridID, 2, 3, 4, 5);

  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, 6);
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, 7);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, "myReference");

  for (i = 0; i < nlon * nlat; i++) area_vec[i] = 0.1 * i;
  gridDefArea(gridID, ap);
  for (i = 0; i < nlon * nlat; i++) mask_vec[i] = i;
  gridDefReducedPoints(gridID, nlon * nlat, mp);
  gridDefComplexPacking(gridID, 1);
  {
    unsigned char uuid[CDI_UUID_SIZE];
    cdiCreateUUID(uuid);
    gridDefUUID(gridID, uuid);
  }

  return gridID;
}

static int
defineZaxis(void)
{
  double vct[3] = { 3.0, 3.3, 3.6 };

  int zaxisID = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID, levs);
  zaxisDefLevel(zaxisID, 2, 8507.3);

  cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, "myName");
  cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, "myLongname");
  cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, "myUnits");

  zaxisDefDatatype(zaxisID, DOUBLE_PRECISION);
  cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, 1);
  zaxisDefVct(zaxisID, 3, vct);
  zaxisDefLbounds(zaxisID, &levs[0]);
  zaxisDefUbounds(zaxisID, &levs[0]);
  zaxisDefWeights(zaxisID, &levs[0]);
  {
    unsigned char uuid[CDI_UUID_SIZE];
    cdiCreateUUID(uuid);
    cdiDefKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);
  }

  return zaxisID;
}

static int
defineTaxis(void)
{
  int taxisID = taxisCreate(TAXIS_ABSOLUTE);

  taxisDefType(taxisID, 0);
  taxisDefVdate(taxisID, 1);
  taxisDefVtime(taxisID, 2);
  taxisDefRdate(taxisID, 3);
  taxisDefRtime(taxisID, 4);
  taxisDefVdateBounds(taxisID, 5, 6);
  taxisDefVtimeBounds(taxisID, 7, 8);
  taxisDefCalendar(taxisID, 1);
  taxisDefTunit(taxisID, 1);
  taxisDefNumavg(taxisID, 1);

  return taxisID;
}

static void
defineStream(int streamID, int vlistID)
{
  streamDefByteorder(streamID, 1);
  streamDefCompType(streamID, 2);
  streamDefCompLevel(streamID, 3);
  streamDefVlist(streamID, vlistID);
}

struct idPair
{
  int id1, id2;
};

static struct idPair
defineVlist(int gridID, int zaxisID, int taxisID)
{
  int vlistID = CDI_UNDEFID;
  int zaxisID2 = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID = vlistCreate();
  int varID1 = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
  int varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARIABLE);
  vlistDefVarName(vlistID, varID1, "varname1");
  {
    int globfac[] = { 23, 42 };
    cdiDefAttInt(vlistID, varID1, "seer's globule factors", CDI_DATATYPE_INT16, 2, globfac);
  }
  vlistDefVarName(vlistID, varID2, "varname2");
  cdiDefAttTxt(vlistID, varID2, "txt demo", 6, "banana");
  vlistDefTaxis(vlistID, taxisID);
  int vlistID2 = vlistCreate();
  vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARIABLE);
  vlistCopy(vlistID2, vlistID);
  return (struct idPair){ vlistID, vlistID2 };
}

static int
defineInstitute(void)
{
  int instID = institutDef(0, 0, "MYINSTITUTE", "myInstitute");
  return instID;
}

static int
defineModel(int instID)
{
  int modelID = modelDef(instID, 0, "resource_copy");
  return modelID;
}

static int destNamespace;

static int
modelRun(MPI_Comm comm)
{
  char *recvBuffer, *sendBuffer;
  int bufferSize;

#ifdef USE_MPI
  const char *fname = "example_resource_copy_mpi.grb";
  cdiPioSerializeSetMPI();
#else
  const char *fname = "example_resource_copy.grb";
#endif

  int gridID = defineGrid();
  int zaxisID = defineZaxis();
  int taxisID = defineTaxis();
  int instID = defineInstitute();
  defineModel(instID);

  struct idPair temp = defineVlist(gridID, zaxisID, taxisID);
  int vlistID = temp.id1;
  int streamID = streamOpenWrite(fname, CDI_FILETYPE_GRB);
  if (streamID < 0) xabort("Could not open file");
  defineStream(streamID, vlistID);
  vlistDestroy(temp.id1);
  vlistDestroy(temp.id2);

  reshPackBufferCreate(&sendBuffer, &bufferSize, &comm);
  recvBuffer = (char *) Malloc((size_t) bufferSize);
#ifdef USE_MPI
  xmpi(MPI_Sendrecv(sendBuffer, bufferSize, MPI_PACKED, 0, 0, recvBuffer, bufferSize, MPI_PACKED, 0, 0, MPI_COMM_SELF,
                    MPI_STATUS_IGNORE));
#else
  memcpy(recvBuffer, sendBuffer, (size_t) bufferSize);
#endif
  namespaceSetActive(destNamespace);
  reshUnpackResources(recvBuffer, bufferSize, &comm, (cdiPostResUpdateHook) 0);
  free(recvBuffer);
  reshPackBufferDestroy(&sendBuffer);

  int differ = reshListCompare(0, 1);

  namespaceSetActive(0);
  streamClose(streamID);
  return differ;
}

int
main(int argc, char *argv[])
{
  int exitCode = 77;
  MPI_Comm commModel;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  commModel = MPI_COMM_WORLD;
#else
  (void) argc;
  (void) argv;
  commModel = 0;
#endif
  destNamespace = namespaceNew();

  exitCode = modelRun(commModel);

#ifdef USE_MPI
  xmpi(MPI_Finalize());
#endif

  return exitCode;
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
