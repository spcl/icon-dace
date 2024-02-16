#ifndef PIO_INTERFACE_
#define PIO_INTERFACE_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <mpi.h>
#include <yaxt.h>

#include "resource_handle.h"
#include "pio_rpc.h"

void cdiPioBufferPartData(int streamID, int varID, int memtype, const void *data, size_t nmiss, Xt_idxlist partDesc);
void cdiPioBufferPartDataGather(int streamID, int varID, int memtype, const void *data, int numBlocks, const int blocklengths[],
                                const int displacements[], size_t nmiss, Xt_idxlist partDesc);

void pioBufferFuncCall(int streamID, struct winHeaderEntry header, const void *data, valPackFunc dataPackFunc);

struct memCpyDataDesc
{
  const void *obj;
  size_t obj_size;
};

void memcpyPackFunc(void *dataDesc, void *buf, int size, int *pos, void *context);

extern float cdiPIOpartInflate_;

void cdiPioStreamWriteVarPart_(int streamID, int varID, int memtype, const void *data, int nmiss, Xt_idxlist partDesc);

void cdiPioClientStreamWinInit(int streamID);
void cdiPioClientStreamWinCreate(int streamID, struct collSpec *cspec);
void cdiPioClientStreamWinDestroy(int streamID);
bool cdiPioClientStreamNeedsFlush(int streamID);
void cdiPioClientStreamWinPost(int streamID);

void cdiPioSetStreamPartDescPreset(int streamID, size_t nVars, const Xt_idxlist partDesc[], const int conversion[]);
struct partDescPreset cdiPioGetStreamPartDescPreset(int streamID);

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
