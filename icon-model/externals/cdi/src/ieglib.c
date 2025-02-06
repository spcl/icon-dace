#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "ieg.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "exse.h"
#include "swap.h"

static int initIegLib = 0;
static int iegDefaultDprec = 0;

// A version string.
#undef LIBVERSION
#define LIBVERSION 1.5.0
#define XSTRING(x) #x
#define STRING(x) XSTRING(x)
static const char ieg_libvers[] = STRING(LIBVERSION);

const char *
iegLibraryVersion(void)
{
  return ieg_libvers;
}

static int IEG_Debug = 0;  // If set to 1, debugging

static void
iegLibInit(void)
{
  const char *envName = "IEG_PRECISION";

  char *envString = getenv(envName);
  if (envString)
    {
      int nrun = (strlen(envString) == 2) ? 1 : 2;
      int pos = 0;
      while (nrun--)
        {
          switch (tolower((int) envString[pos]))
            {
            case 'r':
              {
                switch ((int) envString[pos + 1])
                  {
                  case '4': iegDefaultDprec = EXSE_SINGLE_PRECISION; break;
                  case '8': iegDefaultDprec = EXSE_DOUBLE_PRECISION; break;
                  default: Warning("Invalid digit in %s: %s", envName, envString);
                  }
                break;
              }
            default:
              {
                Warning("Invalid character in %s: %s", envName, envString);
                break;
              }
            }
          pos += 2;
        }
    }

  initIegLib = 1;
}

void
iegDebug(int debug)
{
  if (debug) Message("debug level %d", debug);
  IEG_Debug = debug;
}

static void
iegInit(iegrec_t *iegp)
{
  iegp->checked = 0;
  iegp->byteswap = 0;
  iegp->dprec = 0;
  iegp->refval = 0;
  iegp->datasize = 0;
  iegp->buffersize = 0;
  iegp->buffer = NULL;
}

void
iegInitMem(void *ieg)
{
  iegrec_t *iegp = (iegrec_t *) ieg;
  memset(iegp->ipdb, 0, sizeof(iegp->ipdb));
  memset(iegp->igdb, 0, sizeof(iegp->igdb));
  memset(iegp->vct, 0, sizeof(iegp->vct));
}

void *
iegNew(void)
{
  if (!initIegLib) iegLibInit();

  iegrec_t *iegp = (iegrec_t *) Malloc(sizeof(iegrec_t));
  iegInit(iegp);
  iegInitMem(iegp);

  return (void *) iegp;
}

void
iegDelete(void *ieg)
{
  iegrec_t *iegp = (iegrec_t *) ieg;

  if (iegp)
    {
      if (iegp->buffer) Free(iegp->buffer);
      Free(iegp);
    }
}

int
iegCheckFiletype(int fileID, int *swap)
{
  size_t data = 0;
  size_t dimx = 0, dimy = 0;
  size_t fact = 0;
  unsigned char buffer[1048], *pbuf;

  if (fileRead(fileID, buffer, 4) != 4) return 0;

  size_t blocklen = get_UINT32(buffer);
  size_t sblocklen = get_SUINT32(buffer);

  if (IEG_Debug) Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  // clang-format off
  if (blocklen == 636 || blocklen == 640)
    {
     *swap = 0;
      fact = 4;
      if (fileRead(fileID, buffer, blocklen+8) != blocklen+8) return 0;
      pbuf = buffer+(37+4)*4;    dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+(37+5)*4;    dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if (blocklen == 1040 || blocklen == 1036)
    {
     *swap = 0;
      fact = 8;
      if (fileRead(fileID, buffer, blocklen+8) != blocklen+8) return 0;
      pbuf = buffer+(37+4)*4;    dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+(37+5)*4;    dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if (sblocklen == 636 || sblocklen == 640)
    {
     *swap = 1;
      fact = 4;
      if (fileRead(fileID, buffer, sblocklen+8) != sblocklen+8) return 0;
      pbuf = buffer+(37+4)*4;     dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+(37+5)*4;     dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }
  else if (sblocklen == 1040 || sblocklen == 1036)
    {
     *swap = 1;
      fact = 8;
      if (fileRead(fileID, buffer, sblocklen+8) != sblocklen+8) return 0;
      pbuf = buffer+(37+4)*4;     dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+(37+5)*4;     dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }
  // clang-format on

  fileRewind(fileID);

  if (IEG_Debug) Message("swap = %d fact = %d", *swap, fact);
  if (IEG_Debug) Message("dimx = %lu dimy = %lu data = %lu", dimx, dimy, data);

  int found = data && (dimx * dimy * fact == data || dimx * dimy * 8 == data);
  return found;
}

void
iegCopyMeta(void *dieg, void *sieg)
{
  iegrec_t *diegp = (iegrec_t *) dieg;
  iegrec_t *siegp = (iegrec_t *) sieg;

  // diegp->byteswap = siegp->byteswap;
  diegp->dprec = siegp->dprec;
  diegp->refval = siegp->refval;

  memcpy(diegp->ipdb, siegp->ipdb, sizeof(siegp->ipdb));
  memcpy(diegp->igdb, siegp->igdb, sizeof(siegp->igdb));
  memcpy(diegp->vct, siegp->vct, sizeof(siegp->vct));
}

static int
iegInqData(void *ieg, int prec, void *data)
{
  iegrec_t *iegp = (iegrec_t *) ieg;
  int ierr = 0;
  int byteswap = iegp->byteswap;
  size_t datasize = iegp->datasize;
  void *buffer = iegp->buffer;
  int dprec = iegp->dprec;

  switch (dprec)
    {
    case EXSE_SINGLE_PRECISION:
      {
        if (sizeof(FLT32) == 4)
          {
            if (byteswap) swap4byte(buffer, datasize);

            if (dprec == prec)
              memcpy(data, buffer, datasize * sizeof(FLT32));
            else
              {
                const float *restrict p = (float *) buffer;
                double *restrict q = (double *) data;
                for (size_t i = 0; i < datasize; i++) q[i] = p[i];
              }
          }
        else
          {
            Error("not implemented for %d byte float", sizeof(FLT32));
          }
        break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
        if (sizeof(FLT64) == 8)
          {
            if (byteswap) swap8byte(buffer, datasize);

            if (dprec == prec)
              memcpy(data, buffer, datasize * sizeof(FLT64));
            else
              {
                const double *restrict p = (double *) buffer;
                float *restrict q = (float *) data;
                for (size_t i = 0; i < datasize; i++) q[i] = (float) p[i];
              }
          }
        else
          {
            Error("not implemented for %d byte float", sizeof(FLT64));
          }
        break;
      }
    default:
      {
        Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return ierr;
}

int
iegInqDataSP(void *ieg, float *data)
{
  return iegInqData(ieg, EXSE_SINGLE_PRECISION, (void *) data);
}

int
iegInqDataDP(void *ieg, double *data)
{
  return iegInqData(ieg, EXSE_DOUBLE_PRECISION, (void *) data);
}

static int
iegDefData(iegrec_t *iegp, int prec, const void *data)
{
  int dprec = iegDefaultDprec ? iegDefaultDprec : iegp->dprec;
  iegp->dprec = dprec ? dprec : prec;

  size_t datasize = (size_t) IEG_G_NumLon(iegp->igdb) * (size_t) IEG_G_NumLat(iegp->igdb);
  size_t blocklen = datasize * (size_t) dprec;

  iegp->datasize = datasize;

  if (iegp->buffersize != blocklen)
    {
      iegp->buffersize = blocklen;
      iegp->buffer = Realloc(iegp->buffer, iegp->buffersize);
    }

  switch (dprec)
    {
    case EXSE_SINGLE_PRECISION:
      {
        if (dprec == prec)
          memcpy(iegp->buffer, data, datasize * sizeof(FLT32));
        else
          {
            const double *restrict p = (const double *) data;
            float *restrict q = (float *) iegp->buffer;
            for (size_t i = 0; i < datasize; i++) q[i] = (float) p[i];
          }
        break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
        if (dprec == prec)
          memcpy(iegp->buffer, data, datasize * sizeof(FLT64));
        else
          {
            const float *restrict p = (const float *) data;
            double *restrict q = (double *) iegp->buffer;
            for (size_t i = 0; i < datasize; i++) q[i] = p[i];
          }
        break;
      }
    default:
      {
        Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return 0;
}

int
iegDefDataSP(void *ieg, const float *data)
{
  return iegDefData((iegrec_t *) ieg, EXSE_SINGLE_PRECISION, (void *) data);
}

int
iegDefDataDP(void *ieg, const double *data)
{
  return iegDefData((iegrec_t *) ieg, EXSE_DOUBLE_PRECISION, (void *) data);
}

int
iegRead(int fileID, void *ieg)
{
  iegrec_t *iegp = (iegrec_t *) ieg;
  union
  {
    double d[200];
    float f[200];
    int32_t i32[200];
  } buf;

  if (!iegp->checked)
    {
      int status = iegCheckFiletype(fileID, &iegp->byteswap);
      if (status == 0) Error("Not a IEG file!");
      iegp->checked = 1;
    }

  int byteswap = iegp->byteswap;

  // read header record
  size_t blocklen = binReadF77Block(fileID, byteswap);

  if (fileEOF(fileID)) return -1;

  if (IEG_Debug) Message("blocklen = %lu", blocklen);

  int dprec = 0;
  if (blocklen == 636 || blocklen == 640)
    dprec = 4;
  else if (blocklen == 1040 || blocklen == 1036)
    dprec = 8;
  else
    {
      Warning("unexpecteted header size %d!", (int) blocklen);
      return -1;
    }

  iegp->dprec = dprec;

  binReadInt32(fileID, byteswap, 37, buf.i32);
  for (int i = 0; i < 37; i++) iegp->ipdb[i] = (int) buf.i32[i];

  binReadInt32(fileID, byteswap, 18, buf.i32);
  for (int i = 0; i < 18; i++) iegp->igdb[i] = (int) buf.i32[i];

  if (blocklen == 636 || blocklen == 1036)
    {
      fileRead(fileID, buf.f, 4);
      if (byteswap) swap4byte(buf.f, 1);
      iegp->refval = (double) buf.f[0];
    }
  else
    {
      fileRead(fileID, buf.d, 8);
      if (byteswap) swap8byte(buf.d, 1);
      iegp->refval = (double) buf.d[0];
    }

  binReadInt32(fileID, byteswap, 3, buf.i32);
  for (int i = 0; i < 3; i++) iegp->igdb[18 + i] = (int) buf.i32[i];

  if (dprec == EXSE_SINGLE_PRECISION)
    {
      fileRead(fileID, buf.f, 400);
      if (byteswap) swap4byte(buf.f, 100);
      for (int i = 0; i < 100; i++) iegp->vct[i] = (double) buf.f[i];
    }
  else
    {
      fileRead(fileID, buf.d, 800);
      if (byteswap) swap8byte(buf.d, 100);
      for (int i = 0; i < 100; i++) iegp->vct[i] = buf.d[i];
    }

  size_t blocklen2 = binReadF77Block(fileID, byteswap);

  if (blocklen2 != blocklen)
    {
      Warning("header blocklen differ!");
      return -1;
    }

  iegp->datasize = (size_t) IEG_G_NumLon(iegp->igdb) * (size_t) IEG_G_NumLat(iegp->igdb);

  if (IEG_Debug) Message("datasize = %zu", iegp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  if (iegp->buffersize < blocklen)
    {
      iegp->buffer = Realloc(iegp->buffer, blocklen);
      iegp->buffersize = blocklen;
    }

  if (dprec != (int) (blocklen / iegp->datasize))
    {
      Warning("data precision differ! (h = %d; d = %d)", (int) dprec, (int) (blocklen / iegp->datasize));
      return -1;
    }

  fileRead(fileID, iegp->buffer, blocklen);

  blocklen2 = binReadF77Block(fileID, byteswap);

  if (blocklen2 != blocklen)
    {
      Warning("data blocklen differ!");
      return -1;
    }

  return 0;
}

int
iegWrite(int fileID, void *ieg)
{
  iegrec_t *iegp = (iegrec_t *) ieg;
  union
  {
    INT32 i32[200];
    float fvct[100];
  } buf;
  int dprec = iegp->dprec;
  int byteswap = iegp->byteswap;

  // write header record
  size_t blocklen = (dprec == EXSE_SINGLE_PRECISION) ? 636 : 1040;

  binWriteF77Block(fileID, byteswap, blocklen);

  for (int i = 0; i < 37; i++) buf.i32[i] = (INT32) iegp->ipdb[i];
  binWriteInt32(fileID, byteswap, 37, buf.i32);

  for (int i = 0; i < 18; i++) buf.i32[i] = (INT32) iegp->igdb[i];
  binWriteInt32(fileID, byteswap, 18, buf.i32);

  FLT64 refval = (FLT64) iegp->refval;
  FLT32 refvalf = (FLT32) iegp->refval;
  if (dprec == EXSE_SINGLE_PRECISION)
    binWriteFlt32(fileID, byteswap, 1, &refvalf);
  else
    binWriteFlt64(fileID, byteswap, 1, &refval);

  for (int i = 0; i < 3; i++) buf.i32[i] = (INT32) iegp->igdb[18 + i];
  binWriteInt32(fileID, byteswap, 3, buf.i32);

  if (dprec == EXSE_SINGLE_PRECISION)
    {
      for (int i = 0; i < 100; i++) buf.fvct[i] = (float) iegp->vct[i];
      binWriteFlt32(fileID, byteswap, 100, buf.fvct);
    }
  else
    {
      binWriteFlt64(fileID, byteswap, 100, iegp->vct);
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  iegp->datasize = (size_t) iegp->igdb[4] * (size_t) iegp->igdb[5];
  blocklen = iegp->datasize * (size_t) dprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch (dprec)
    {
    case EXSE_SINGLE_PRECISION:
      {
        binWriteFlt32(fileID, byteswap, iegp->datasize, (FLT32 *) iegp->buffer);
        break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
        binWriteFlt64(fileID, byteswap, iegp->datasize, (FLT64 *) iegp->buffer);
        break;
      }
    default:
      {
        Error("unexpected data precision %d", dprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

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
