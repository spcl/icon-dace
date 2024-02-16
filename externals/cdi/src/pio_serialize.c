#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>

#include <mpi.h>

#include "cdi.h"
#include "namespace.h"
#include "pio_serialize.h"
#include "pio_util.h"

#if MPI_VERSION == 1 || (MPI_VERSION == 2 && MPI_SUBVERSION < 2)
#define CDI_DT_MATCH_NEEDED 1
#endif

static
#ifndef CDI_DT_MATCH_NEEDED
    const
#endif
    struct
{
  int cdidt;
  MPI_Datatype mpidt;
} dtDict[] = {
#ifdef CDI_DT_MATCH_NEEDED
  { CDI_DATATYPE_INT8, MPI_SIGNED_CHAR }, { CDI_DATATYPE_INT16, MPI_SHORT },   { CDI_DATATYPE_INT32, MPI_INT },
  { CDI_DATATYPE_UINT32, MPI_INT },
#else
  { CDI_DATATYPE_INT8, MPI_INT8_T },     { CDI_DATATYPE_INT16, MPI_INT16_T }, { CDI_DATATYPE_INT32, MPI_INT32_T },
  { CDI_DATATYPE_UINT32, MPI_UINT32_T },
#endif
  { CDI_DATATYPE_INT, MPI_INT },          { CDI_DATATYPE_UINT, MPI_UNSIGNED }, { CDI_DATATYPE_FLT64, MPI_DOUBLE },
  { CDI_DATATYPE_FLT, MPI_DOUBLE },       { CDI_DATATYPE_TXT, MPI_CHAR },      { CDI_DATATYPE_UCHAR, MPI_UNSIGNED_CHAR },
  { CDI_DATATYPE_LONG, MPI_LONG },
};

static inline size_t
lookupDt(int datatype)
{
  for (size_t i = 0; i < sizeof(dtDict) / sizeof(dtDict[0]); ++i)
    if (dtDict[i].cdidt == datatype) return i;
  abort();
}

#ifdef CDI_DT_MATCH_NEEDED
static int dtDictMatchComplete = 0;

static inline void
dtDictFixMPIType(size_t i, int typeclass, int size)
{
  MPI_Aint lb, extent;
  xmpi(MPI_Type_get_extent(dtDict[i].mpidt, &lb, &extent));
  if ((int) extent != size)
    {
      /* type size mismatch needs to be fixed */
      MPI_Type_match_size(typeclass, size, &dtDict[i].mpidt);
    }
}

static void
setupDtDict()
{
  dtDictFixMPIType(lookupDt(CDI_DATATYPE_INT8), MPI_TYPECLASS_INTEGER, (int) sizeof(int8_t));
  dtDictFixMPIType(lookupDt(CDI_DATATYPE_INT16), MPI_TYPECLASS_INTEGER, (int) sizeof(int16_t));
  dtDictFixMPIType(lookupDt(CDI_DATATYPE_INT32), MPI_TYPECLASS_INTEGER, (int) sizeof(int32_t));
  dtDictFixMPIType(lookupDt(CDI_DATATYPE_UINT32), MPI_TYPECLASS_INTEGER, (int) sizeof(uint32_t));
  dtDictMatchComplete = 1;
}
#endif

static int
serializeGetSizeMPI(int count, int datatype, void *context)
{
  int size;
  xmpi(MPI_Pack_size(count, dtDict[lookupDt(datatype)].mpidt, *(MPI_Comm *) context, &size));
  return size;
}

static void
serializePackMPI(void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context)
{
  xmpi(MPI_Pack(data, count, dtDict[lookupDt(datatype)].mpidt, buf, buf_size, position, *(MPI_Comm *) context));
}

static void
serializeUnpackMPI(void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context)
{
  xmpi(MPI_Unpack(buf, buf_size, position, data, count, dtDict[lookupDt(datatype)].mpidt, *(MPI_Comm *) context));
}

void
cdiPioSerializeSetMPI(void)
{
#ifdef CDI_DT_MATCH_NEEDED
  if (!dtDictMatchComplete) setupDtDict();
#endif
  namespaceSwitchSet(NSSWITCH_SERIALIZE_GET_SIZE, NSSW_FUNC(serializeGetSizeMPI));
  namespaceSwitchSet(NSSWITCH_SERIALIZE_PACK, NSSW_FUNC(serializePackMPI));
  namespaceSwitchSet(NSSWITCH_SERIALIZE_UNPACK, NSSW_FUNC(serializeUnpackMPI));
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
