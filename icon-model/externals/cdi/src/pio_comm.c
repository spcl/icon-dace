#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <yaxt.h>

#include "cdi.h"
#include "dmemory.h"

#include "pio.h"
#include "pio_comm.h"
#include "pio_util.h"

#include "cdipio.h"

struct cdiPioComms
{
  int IOMode;
  int isProcIO;

  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int root;

  MPI_Comm commPio; /* intracommunicator of PIO servers */
  int sizePio;
  int rankPio;

  MPI_Comm commColl; /* intracommunicator of collector
                      * subset of PIO servers only */
  int sizeColl;
  int rankColl;

  MPI_Comm commCalc;
  int rankCalc, sizeCalc;

  /**
   * communicator connecting group of clients and group of collectors
   * notice: it is always implicit that the group of collectors is no
   * larger than the group of clients
   */
  MPI_Comm pioInterComm;
  /**
   * group containing only the active remote
   * processes, i.e.
   * - on collector processes contains all processes that are a client of
   *   the server
   * - on client processes contains only the server process used for
   *   data/RPC communication with this client process
   */
  MPI_Group remoteGroup;
  /**
   * Communicator containing both clients and collectors, needed since
   * MPI_Win_create only accepts intracommunicators.
   * Constructed such that collectors occupy the top ranks in the same
   * sequence they have on their own communicator.
   */
  MPI_Comm collClientIntraComm;
};

static void setupCollectorInterComm(struct cdiPioComms *p, int isCollector, MPI_Comm localGroupComm, int globRankOfGroupLeader[2]);

static void
pioInfoInit(struct cdiPioComms *p, MPI_Comm commGlob, int clientServerRole)
{
  xassert(p && commGlob != MPI_COMM_NULL);

  xmpi(MPI_Comm_dup(commGlob, &p->commGlob));
  xmpi(MPI_Comm_size(p->commGlob, &p->sizeGlob));
  xmpi(MPI_Comm_rank(p->commGlob, &p->rankGlob));
  p->root = 0;

  p->commColl = MPI_COMM_NULL;
  p->rankColl = -1;

  int isProcIO = (clientServerRole != PIO_ROLE_CLIENT), isCollector = (clientServerRole & PIO_ROLE_COLLECTOR);
  p->isProcIO = isProcIO;

  {
    MPI_Comm temp;
    /* key is 0 on all clients
     *        0 on all servers that are also collectors
     *        1 on all servers that are not collectors
     */
    int key = isProcIO - isCollector;
    xmpi(MPI_Comm_split(p->commGlob, isProcIO, key, &temp));
    enum
    {
      PIO = 0,
      Client = 1,
      numGroups
    };
    int globRankOfGroupLeader[numGroups] = { 0, 0 };
    if (isProcIO)
      {
        p->commPio = temp;
        xmpi(MPI_Comm_rank(temp, &p->rankPio));
        p->commCalc = MPI_COMM_NULL;
        p->rankCalc = -1;
        if (p->rankPio == 0) globRankOfGroupLeader[PIO] = p->rankGlob;
      }
    else
      {
        p->commPio = MPI_COMM_NULL;
        p->rankPio = -1;
        p->commCalc = temp;
        xmpi(MPI_Comm_rank(temp, &p->rankCalc));
        xmpi(MPI_Comm_size(temp, &p->sizeCalc));
        if (p->rankCalc == 0) globRankOfGroupLeader[Client] = p->rankGlob;
      }
    int nProcsIO;
    {
      enum
      {
        numGlobSums = 4
      };
      int iBuf[numGlobSums]
          = { isProcIO, (clientServerRole & PIO_ROLE_COLLECTOR), globRankOfGroupLeader[0], globRankOfGroupLeader[1] };
      xmpi(MPI_Allreduce(MPI_IN_PLACE, iBuf, numGlobSums, MPI_INT, MPI_SUM, commGlob));
      p->sizePio = nProcsIO = iBuf[0];
      p->sizeCalc = p->sizeGlob - nProcsIO;
      p->sizeColl = iBuf[1];
      globRankOfGroupLeader[0] = iBuf[2];
      globRankOfGroupLeader[1] = iBuf[3];
    }
    if (isProcIO)
      {
        MPI_Group pioGroup, collGroup;
        int collRange[1][3] = { { 0, p->sizeColl - 1, 1 } };
        xmpi(MPI_Comm_group(p->commPio, &pioGroup));
        xmpi(MPI_Group_range_incl(pioGroup, 1, collRange, &collGroup));
        xmpi(MPI_Comm_create(p->commPio, collGroup, &temp));
        if (isCollector)
          {
            xmpi(MPI_Comm_rank(temp, &p->rankColl));
            xt_mpi_comm_mark_exclusive(temp);
          }
        xmpi(MPI_Group_free(&pioGroup));
        xmpi(MPI_Group_free(&collGroup));
        p->commColl = temp;
      }
    if (isCollector || !isProcIO)
      setupCollectorInterComm(p, isCollector, temp, globRankOfGroupLeader);
    else /* isProcIO && !isCollector */
      {
        p->pioInterComm = MPI_COMM_NULL;
        p->collClientIntraComm = MPI_COMM_NULL;
        p->remoteGroup = MPI_GROUP_NULL;
      }

    /* match output method and process distribution */
    int IOMode = p->IOMode, sizeGlob = p->sizeGlob;
    if (((IOMode != PIO_NONE && (nProcsIO <= 0 || nProcsIO > sizeGlob - 1))) || (IOMode == PIO_NONE && nProcsIO != 1))
      xabort("DISTRIBUTION OF TASKS ON PROCS IS NOT VALID.\n"
             "nProcsIO=%d, sizeGlob=%d\n",
             nProcsIO, sizeGlob);
  }
}

/* create client <-> collector intercomm and fill-in related information */
static void
setupCollectorInterComm(struct cdiPioComms *p, int isCollector, MPI_Comm localGroupComm, int globRankOfGroupLeader[2])
{
  xmpi(MPI_Intercomm_create(localGroupComm, 0, p->commGlob, globRankOfGroupLeader[isCollector], 1, &p->pioInterComm));
  xmpi(MPI_Intercomm_merge(p->pioInterComm, isCollector, &p->collClientIntraComm));
  MPI_Group collClientGroup;
  int remoteRanksRange[3];
  xmpi(MPI_Comm_group(p->collClientIntraComm, &collClientGroup));
  int numClients = p->sizeCalc, numColl = p->sizeColl;
  if (isCollector)
    {
      int collRank = p->rankColl;
      remoteRanksRange[0] = cdiPioClientRangeStart(collRank, numClients, numColl);
      remoteRanksRange[1] = cdiPioClientRangeStart(collRank + 1, numClients, numColl) - 1;
    }
  else /* !isProcIO i.e. is client */
    {
      int clientRank = p->rankCalc;
      int collRank = cdiPioCollRank(clientRank, numClients, numColl);
      remoteRanksRange[1] = remoteRanksRange[0] = numClients + collRank;
    }
  remoteRanksRange[2] = 1;
  xmpi(MPI_Group_range_incl(collClientGroup, 1, &remoteRanksRange, &p->remoteGroup));
  xmpi(MPI_Group_free(&collClientGroup));
}

void
cdiPioCommInit(MPI_Comm commGlob, int IOMode, int clientServerRole)
{
  xassert(IOMode >= PIO_MINIOMODE && IOMode <= PIO_MAXIOMODE);
  if (cdiPioExtraNSKeys[cdiPioEKComms] == 0) cdiPioExtraNSKeys[cdiPioEKComms] = cdiNamespaceSwitchNewKey();
  struct cdiPioComms *info = Malloc(sizeof(*info));
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKComms], NSSW_DATA(info));
  info->IOMode = IOMode;
  pioInfoInit(info, commGlob, clientServerRole);
}

void
cdiPioCommFinalize(void)
{
  xassert(cdiPioExtraNSKeys[cdiPioEKComms] != 0);
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;

  {
    MPI_Comm *comms[]
        = { &info->commGlob, &info->commPio, &info->commColl, &info->commCalc, &info->pioInterComm, &info->collClientIntraComm };

    for (size_t i = 0; i < sizeof(comms) / sizeof(comms[0]); ++i)
      if (*comms[i] != MPI_COMM_NULL) xmpi(MPI_Comm_free(comms[i]));
  }

  {
    MPI_Group *groups[] = { &info->remoteGroup };
    for (size_t i = 0; i < sizeof(groups) / sizeof(groups[0]); ++i)
      if (*groups[i] != MPI_GROUP_NULL) xmpi(MPI_Group_free(groups[i]));
  }

  Free(info);
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKComms], NSSW_DATA(NULL));
}

MPI_Comm
commInqCommGlob(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->commGlob != MPI_COMM_NULL);
  return info->commGlob;
}

int
commInqSizeGlob(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->sizeGlob != CDI_UNDEFID);
  return info->sizeGlob;
}

int
commInqRankGlob(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->rankGlob != -1);
  return info->rankGlob;
}

MPI_Group
cdiPioInqRemoteGroup(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info);
  return info->remoteGroup;
}

int
commInqIsProcIO(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->isProcIO != CDI_UNDEFID);
  return info->isProcIO;
}

int
commInqIOMode(void)
{
  struct cdiPioComms *info
      = cdiPioExtraNSKeys[cdiPioEKComms] == 0 ? NULL : namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  return info != NULL ? info->IOMode : PIO_NONE;
}

MPI_Comm
commInqCommPio(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->commPio != MPI_COMM_NULL && info->isProcIO == 1);
  return info->commPio;
}

int
commInqSizePio(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->sizePio > 0);
  return info->sizePio;
}

MPI_Comm
commInqCommModel(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->commCalc != MPI_COMM_NULL && info->isProcIO == 0);
  return info->commCalc;
}

int
commInqRankPio(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->rankPio >= 0 && info->isProcIO == 1);
  return info->rankPio;
}

int
commInqRankModel(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->rankCalc >= 0 && info->isProcIO == 0);
  return info->rankCalc;
}

int
cdiPioCommInqSizeClients(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->sizeCalc >= 0);
  return info->sizeCalc;
}

MPI_Comm
commInqCommColl(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->commColl != MPI_COMM_NULL);
  return info->commColl;
}

int
commInqSizeColl(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->sizeColl >= 0);
  return info->sizeColl;
}

int
commInqRankColl(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info);
  return info->rankColl;
}

MPI_Comm
cdiPioInqInterComm(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info);
  return info->pioInterComm;
}

MPI_Comm
cdiPioInqCollClientIntraComm(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info);
  return info->collClientIntraComm;
}

int
commInqSpecialRank(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info);
  return info->sizeColl;
}

MPI_Comm
commInqCommCalc(void)
{
  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;
  xassert(info && info->commCalc != MPI_COMM_NULL);
  return info->commCalc;
}

void
commPrint(FILE *fp)
{
  if (ddebug == 0) return;

  struct cdiPioComms *info = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKComms]).data;

  xassert(info != NULL);

  fprintf(fp, "\n");
  fprintf(fp, "######## pioinfo PE%d ###########\n", info->rankGlob);
  fprintf(fp, "#\n");
  fprintf(fp, "# IOMode      = %d\n", info->IOMode);
  fprintf(fp, "# isProcIO    = %d\n", info->isProcIO);
  fprintf(fp, "#\n");
  fprintf(fp, "# commGlob    = %d\n", (int) MPI_Comm_c2f(info->commGlob));
  fprintf(fp, "# sizeGlob    = %d\n", info->sizeGlob);
  fprintf(fp, "# rankGlob    = %d\n", info->rankGlob);
  fprintf(fp, "#\n");
  fprintf(fp, "# commPio     = %d\n", (int) MPI_Comm_c2f(info->commPio));
  fprintf(fp, "# sizePio     = %d\n", info->sizePio);
  fprintf(fp, "# rankPio     = %d\n", info->rankPio);
  fprintf(fp, "#\n");
  fprintf(fp, "# commColl    = %d\n", (int) MPI_Comm_c2f(info->commColl));
  fprintf(fp, "# sizeColl    = %d\n", info->sizeColl);
  fprintf(fp, "# rankColl    = %d\n", info->rankColl);
  fprintf(fp, "#\n");
  fprintf(fp, "# commCalc    = %d\n", (int) MPI_Comm_c2f(info->commCalc));
  fprintf(fp, "# sizeCalc  = %d\n", info->sizeCalc);
  fprintf(fp, "# rankCalc  = %d\n", info->rankCalc);
  fprintf(fp, "#\n");
  fprintf(fp, "############################\n");
  fprintf(fp, "\n");
}

/************************************************************************/
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
