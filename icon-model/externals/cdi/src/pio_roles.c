#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include "cdipio.h"
#include "error.h"
#include "pio_util.h"

static int cdiPioCSRSelectRange(MPI_Comm commSuper, int IOMode, int clientRangeStart, int clientRangeEnd, int specialRank);

int
cdiPioCSRLastN(MPI_Comm commSuper, int IOMode, int nProcsIO)
{
  int commSize;
  xmpi(MPI_Comm_size(commSuper, &commSize));
  return cdiPioCSRSelectRange(commSuper, IOMode, 0, commSize - nProcsIO - 1, commSize - 1);
}

int
cdiPioCSRFirstN(MPI_Comm commSuper, int IOMode, int nProcsIO)
{
  int commSize;
  xmpi(MPI_Comm_size(commSuper, &commSize));
  return cdiPioCSRSelectRange(commSuper, IOMode, nProcsIO, commSize - 1, nProcsIO - 1);
}

static int
cdiPioCSRSelectRange(MPI_Comm commSuper, int IOMode, int clientRangeStart, int clientRangeEnd, int specialRank)
{
  int commRank, commSize;
  xmpi(MPI_Comm_size(commSuper, &commSize));
  xmpi(MPI_Comm_rank(commSuper, &commRank));
  int role;
  switch (IOMode)
    {
    case PIO_NONE:
      xassert(commRank == 0 && commSize == 1);
      role = PIO_NONE;
      break;
    case PIO_MPI:
    case PIO_MPI_FW_ORDERED:
    case PIO_MPI_FW_AT_ALL:
    case PIO_MPI_FW_AT_REBLOCK:
      role = (commRank >= clientRangeStart && commRank <= clientRangeEnd) ? PIO_ROLE_CLIENT : PIO_ROLE_WRITER_COLLECTOR;
      break;
    case PIO_WRITER:
    case PIO_ASYNCH:
      role = (commRank >= clientRangeStart && commRank <= clientRangeEnd)
                 ? PIO_ROLE_CLIENT
                 : ((commRank == specialRank) ? PIO_ROLE_WRITER : PIO_ROLE_COLLECTOR);
      break;
    case PIO_FPGUARD:
      role = (commRank >= clientRangeStart && commRank <= clientRangeEnd)
                 ? PIO_ROLE_CLIENT
                 : ((commRank == specialRank) ? PIO_ROLE_FPGUARD : PIO_ROLE_WRITER_COLLECTOR);
      break;
    default: xabort("Invalid mode requested %d", IOMode);
    }
  return role;
}

int
cdiPioCSRBalanced(MPI_Comm commSuper, int IOMode, int nProcsIO)
{
  int commRank, commSize, numClients, clientsPerCollectorMin, clientsPerCollectorMax, rest;
  xmpi(MPI_Comm_size(commSuper, &commSize));
  xmpi(MPI_Comm_rank(commSuper, &commRank));
  numClients = commSize - nProcsIO;
  int role, specialRole = PIO_ROLE_FPGUARD, collType = PIO_ROLE_WRITER_COLLECTOR;
  switch (IOMode)
    {
    case PIO_NONE:
      xassert(commRank == 0 && commSize == 1);
      role = PIO_NONE;
      break;
    case PIO_MPI:
    case PIO_MPI_FW_ORDERED:
    case PIO_MPI_FW_AT_ALL:
    case PIO_MPI_FW_AT_REBLOCK:
      clientsPerCollectorMin = numClients / nProcsIO;
      rest = numClients % nProcsIO;
      clientsPerCollectorMax = clientsPerCollectorMin + (rest != 0);
      if (commRank == commSize - 1)
        role = PIO_ROLE_WRITER_COLLECTOR;
      else if (commRank <= (clientsPerCollectorMax + 1) * (nProcsIO - rest))
        role = ((commRank + 1) % (clientsPerCollectorMax + 1)) != 0 ? PIO_ROLE_CLIENT : PIO_ROLE_WRITER_COLLECTOR;
      else
        role = ((commRank - clientsPerCollectorMax * (nProcsIO - rest) + 1) % (clientsPerCollectorMin + 1)) != 0
                   ? PIO_ROLE_CLIENT
                   : PIO_ROLE_WRITER_COLLECTOR;
      break;
    case PIO_WRITER:
    case PIO_ASYNCH:
      specialRole = PIO_ROLE_WRITER;
      collType = PIO_ROLE_COLLECTOR;
      /*-fallthrough*/
    case PIO_FPGUARD:
      clientsPerCollectorMin = numClients / (nProcsIO - 1);
      rest = numClients % nProcsIO;
      clientsPerCollectorMax = clientsPerCollectorMin + (rest != 0);
      if (commRank == commSize - 1)
        role = specialRole;
      else if (commRank <= clientsPerCollectorMax * (nProcsIO - 1 - rest))
        role = ((commRank + 1) % (clientsPerCollectorMax + 1)) != 0 ? PIO_ROLE_CLIENT : collType;
      else
        role = ((commRank - clientsPerCollectorMax * (nProcsIO - 1 - rest) + 1) % (clientsPerCollectorMin + 1)) != 0
                   ? PIO_ROLE_CLIENT
                   : collType;
      break;
    default: xabort("Invalid mode requested %d", IOMode);
    }
  return role;
}
