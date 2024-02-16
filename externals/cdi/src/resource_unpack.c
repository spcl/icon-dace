#include "cdi.h"
#include "dmemory.h"
#include "grid.h"
#include "institution.h"
#include "model.h"
#include "cdi_int.h"
#include "vlist.h"
#include "namespace.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "taxis.h"
#include "zaxis.h"

int (*reshDistGridUnpack)(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                          int force_id);

/*****************************************************************************/

int
reshUnpackResources(char *unpackBuffer, int unpackBufferSize, void *context, cdiPostResUpdateHook postHook)
{
  int updateType, resH, originNamespace;
  int unpackBufferPos = 0;
  int numAssociations = 0, sizeAssociations = 16;
  struct streamAssoc *associations = (struct streamAssoc *) Malloc(sizeof(associations[0]) * (size_t) sizeAssociations);

  {
    int msgHdr[2];
    serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, &msgHdr, 2, CDI_DATATYPE_INT, context);
    if (msgHdr[0] != START) xabort("error parsing resource serialization buffer");
    originNamespace = msgHdr[1];
  }
  while (unpackBufferPos < unpackBufferSize)
    {
      serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, &updateType, 1, CDI_DATATYPE_INT, context);
      if (updateType == END) break;
      int updatedResH;
      switch (updateType)
        {
        case GRID: updatedResH = gridUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1); break;
        case ZAXIS: updatedResH = zaxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1); break;
        case TAXIS: updatedResH = taxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1); break;
        case INSTITUTE:
          updatedResH = instituteUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1);
          break;
        case MODEL: updatedResH = modelUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1); break;
        case STREAM:
          if (sizeAssociations == numAssociations)
            associations = (struct streamAssoc *) Realloc(associations, sizeof(associations[0]) * (size_t) (sizeAssociations *= 2));
          {
            struct streamAssoc newAssoc;
            associations[numAssociations] = newAssoc
                = streamUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context);
            updatedResH = newAssoc.streamID;
          }
          ++numAssociations;
          break;
        case VLIST: updatedResH = vlistUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1); break;
        case DIST_GRID:
          updatedResH = reshDistGridUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, originNamespace, context, 1);
          break;
        case RESH_DELETE:
          serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos, &resH, 1, CDI_DATATYPE_INT, context);
          resH = namespaceAdaptKey(resH, originNamespace);
          reshDestroy(resH);
          updatedResH = resH;
          reshSetStatus(resH, NULL, RESH_UNUSED);
          break;
        default: xabort("Invalid/unexpected serialization type %d or transfer error!", updateType);
        }
      if (postHook != (cdiPostResUpdateHook) 0) postHook(updatedResH, updateType);
    }
  for (int i = 0; i < numAssociations; ++i)
    {
      cdiStreamSetupVlist(stream_to_pointer(associations[i].streamID), namespaceAdaptKey(associations[i].vlistID, originNamespace));
    }
  Free(associations);
  return unpackBufferPos;
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
