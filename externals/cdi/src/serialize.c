#include <inttypes.h>
#include <limits.h>
#include <string.h>

#include "cdi.h"
#include "error.h"
#include "serialize.h"
#include "namespace.h"

int
serializeGetSize(int count, int datatype, void *context)
{
  int (*serialize_get_size_p)(int count, int datatype, void *context)
      = (int (*)(int, int, void *)) namespaceSwitchGet(NSSWITCH_SERIALIZE_GET_SIZE).func;
  return serialize_get_size_p(count, datatype, context);
}

void
serializePack(const void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context)
{
  void (*serialize_pack_p)(const void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context)
      = (void (*)(const void *, int, int, void *, int, int *, void *)) namespaceSwitchGet(NSSWITCH_SERIALIZE_PACK).func;
  serialize_pack_p(data, count, datatype, buf, buf_size, position, context);
}

void
serializeUnpack(const void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context)
{
  void (*serialize_unpack_p)(const void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context)
      = (void (*)(const void *, int, int *, void *, int, int, void *)) namespaceSwitchGet(NSSWITCH_SERIALIZE_UNPACK).func;
  serialize_unpack_p(buf, buf_size, position, data, count, datatype, context);
}

int
serializeGetSizeInCore(int count, int datatype, void *context)
{
  int elemSize;
  (void) context;
  switch (datatype)
    {
    case CDI_DATATYPE_INT8: elemSize = sizeof(int8_t); break;
    case CDI_DATATYPE_INT16: elemSize = sizeof(int16_t); break;
    case CDI_DATATYPE_UINT32: elemSize = sizeof(uint32_t); break;
    case CDI_DATATYPE_INT: elemSize = sizeof(int); break;
    case CDI_DATATYPE_UINT: elemSize = sizeof(unsigned); break;
    case CDI_DATATYPE_FLT:
    case CDI_DATATYPE_FLT64: elemSize = sizeof(double); break;
    case CDI_DATATYPE_TXT:
    case CDI_DATATYPE_UCHAR: elemSize = 1; break;
    case CDI_DATATYPE_LONG: elemSize = sizeof(long); break;
    default: xabort("Unexpected datatype");
    }
  return count * elemSize;
}

void
serializePackInCore(const void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context)
{
  int size = serializeGetSize(count, datatype, context);
  int pos = *position;
  xassert(INT_MAX - pos >= size && buf_size - pos >= size);
  memcpy((unsigned char *) buf + pos, data, (size_t) size);
  pos += size;
  *position = pos;
}

void
serializeUnpackInCore(const void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context)
{
  int size = serializeGetSize(count, datatype, context);
  int pos = *position;
  xassert(INT_MAX - pos >= size && buf_size - pos >= size);
  memcpy(data, (unsigned char *) buf + pos, (size_t) size);
  pos += size;
  *position = pos;
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
