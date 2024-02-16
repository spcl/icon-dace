#ifndef PIO_INFO_
#define PIO_INFO_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <mpi.h>
#include <stdio.h>

void cdiPioCommInit(MPI_Comm commGlob, int IOMode, int clientServerRole);
void cdiPioCommFinalize(void);

MPI_Comm commInqCommGlob(void);
int commInqSizeGlob(void);
int commInqRankGlob(void);

/* returns the inter-communicator used to connect from
 * client to server group */
MPI_Comm cdiPioInqInterComm(void);

int commInqIsProcIO(void);
int commInqIOMode(void);

void commDefCommPio(void);
MPI_Comm commInqCommPio(void);
int commInqSizePio(void);
int commInqRankPio(void);

MPI_Comm commInqCommModel(void);
int commInqRankModel(void);
int cdiPioCommInqSizeClients(void);

int commInqSpecialRank(void);

/* query communicator of processes active in collection of data */
MPI_Comm commInqCommColl(void);
int commRankGlob2CollID(int rankGlob);
int commInqSizeColl(void);
int commInqRankColl(void);
MPI_Comm commInqCommCalc(void);

void commPrint(FILE *);

/* query communicator containing clients and collectors,
 * the collectors occupy the top-ranks
 */
MPI_Comm cdiPioInqCollClientIntraComm(void);

/* query group of corresponding processes in other group on
   intra-communicator */
MPI_Group cdiPioInqRemoteGroup(void);

static inline int
cdiPioCollRank(int clientRank, int numClients, int numColl)
{
  return (int) (((long long) clientRank * (long long) numColl + (long long) (numColl - 1)) / numClients);
}

static inline int
cdiPioClientRangeStart(int collRank, int numClients, int numColl)
{
  return (int) (((long long) collRank * (long long) numClients) / (long long) numColl);
}

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
