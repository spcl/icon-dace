#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>

#include <mpi.h>
#include <yaxt.h>
#include "cdi.h"
#include "cdipio.h"
#include "dmemory.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "resource_unpack.h"

extern int reshListCompare(int, int);

enum
{
  IOMode = PIO_NONE,
  nProcsIO = 1,
  DOUBLE_PRECISION = 8,
  nlon = 12,
  nlat = 6,
  nlev = 5,
  ntsteps = 3
};

static double lons[nlon] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
static double lats[nlat] = { -75, -45, -15, 15, 45, 75 };
static double levs[nlev] = { 101300, 92500, 85000, 50000, 20000 };

static int
defineGrid()
{
  int gridID = CDI_UNDEFID;
  int mask_vec[nlon * nlat];
  const int *mp = &mask_vec[0];
  double area_vec[nlon * nlat];
  const double *ap = &area_vec[0];
  int i;

  gridID = gridCreate(GRID_LONLAT, nlon * nlat);
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

  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, 6);
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, 7);
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, "myReference");

  /* gridDefLCC ( gridID, double originLon, double originLat,  */
  /* 	   double lonParY, double lat1, double lat2, double xinc, double yinc, int projflag, int scanflag); */
  /* gridDefLcc2 ( gridID, double earth_radius, double lon_0,  */
  /* 	    double lat_0, double lat_1,double lat_2);*/
  /* gridDefLaea ( gridID, double earth_radius, double lon_0, double lat_0); */
  for (i = 0; i < nlon * nlat; i++) area_vec[i] = 0.1 * i;
  gridDefArea(gridID, ap);
  for (i = 0; i < nlon * nlat; i++) mask_vec[i] = i;
  gridDefReducedPoints(gridID, nlon * nlat, mp);
  gridDefComplexPacking(gridID, 1);

  return gridID;
}

static int
defineZaxis()
{
  int zaxisID = CDI_UNDEFID;
  double vct[3] = { 3.0, 3.3, 3.6 };

  zaxisID = zaxisCreate(ZAXIS_PRESSURE, nlev);
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

  return zaxisID;
}

static int
defineTaxis()
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

static int
defineVlist(int gridID, int zaxisID, int taxisID)
{
  int vlistID = CDI_UNDEFID;
  int zaxisID2 = zaxisCreate(ZAXIS_SURFACE, 1);
  int varID1, varID2;

  vlistID = vlistCreate();
  varID1 = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARIABLE);
  vlistDefVarName(vlistID, varID1, "varname1");
  {
    int globfac[] = { 23, 42 };
    cdiDefAttInt(vlistID, varID1, "seer's globule factors", CDI_DATATYPE_INT16, 2, globfac);
  }
  vlistDefVarName(vlistID, varID2, "varname2");
  cdiDefAttTxt(vlistID, varID2, "txt demo", 6, "banana");
  vlistDefTaxis(vlistID, taxisID);
  return vlistID;
}

static int
defineInstitute()
{
  int instID = CDI_UNDEFID;

  instID = institutDef(0, 0, "MYINSTITUTE", "myInstitute");

  return instID;
}

static void
defineModel(int instID)
{
  modelDef(instID, 0, "myModel");
}

static void
printResources()
{
  FILE *fp = fopen("reshArrayModel.txt", "w");
  if (!fp) xabort("%s", "could not open file");
  reshListPrint(fp);
  fclose(fp);
}

static void
modelRun(MPI_Comm comm)
{
  int gridID, zaxisID, taxisID, instID, vlistID, streamID;

  char *recvBuffer, *sendBuffer;
  int bufferSize, differ;
  MPI_Status status;

  namespaceSetActive(0);

  gridID = defineGrid();
  zaxisID = defineZaxis();
  taxisID = defineTaxis();
  instID = defineInstitute();
  defineModel(instID);
  vlistID = defineVlist(gridID, zaxisID, taxisID);
  streamID = streamOpenWrite("example.grb", CDI_FILETYPE_GRB);
  if (streamID < 0) xabort("Could not open file");
  defineStream(streamID, vlistID);

  reshPackBufferCreate(&sendBuffer, &bufferSize, &comm);
  xmpi(MPI_Send(sendBuffer, bufferSize, MPI_PACKED, 0, 0, comm));
  recvBuffer = Malloc((size_t) bufferSize);
  xmpi(MPI_Recv(recvBuffer, bufferSize, MPI_PACKED, 0, 0, comm, &status));

  namespaceSetActive(1);
  reshUnpackResources(recvBuffer, bufferSize, &comm, (cdiPostResUpdateHook) 0);
  Free(recvBuffer);
  reshPackBufferDestroy(&sendBuffer);

  differ = reshListCompare(0, 1);
  printf("The resource arrays %s.\n", differ ? "differ" : "are equal");
  printResources();

  namespaceSetActive(0);
  streamClose(streamID);
  return;
}

int
main(int argc, char *argv[])
{
  int sizeGlob, pioNamespace;
  MPI_Comm commGlob, commModel;

  MPI_Init(&argc, &argv);
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi(MPI_Comm_set_errhandler(commGlob, MPI_ERRORS_RETURN));
  xmpi(MPI_Comm_size(commGlob, &sizeGlob));

  if (sizeGlob != 1) xabort("test transition of resource array only with 1 PE.");

  if (nProcsIO != 1) xabort("bad distribution of tasks on PEs");

  commModel = pioInit(commGlob, nProcsIO, IOMode, &pioNamespace, 1.0f, cdiPioNoPostCommSetup);
  namespaceSetActive(pioNamespace);

  modelRun(commModel);

  xt_finalize();
  MPI_Finalize();

  return 0;
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
