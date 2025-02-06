#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#ifdef HAVE_PPM_CORE
#include <ppm/ppm.h>
#endif
#endif

#include "cdi.h"
#include "error.h"
#include "dmemory.h"
#include "pio_write.h"
#include "simple_model_helper.h"
#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#endif

static const char default_fname_prefix[] = "example";

static const struct model_config default_setup
#ifdef __cplusplus
    = { 12,
        6,
        3,
        5,
        5,
        CDI_FILETYPE_GRB,
        CDI_DATATYPE_PACK24,
        PIO_WRITE_CONFIG_CHECKSUM_FLAG,
        TAXIS_ABSOLUTE,
        -1,
        "grb",
        default_fname_prefix };
#else
    = {
        .nlon = 12,
        .nlat = 6,
        .nts = 3,
        .max_nlev = 5,
        .nvars = 5,
        .filetype = CDI_FILETYPE_GRB,
        .datatype = CDI_DATATYPE_PACK24,
        .flags = PIO_WRITE_CONFIG_CHECKSUM_FLAG,
        .taxistype = TAXIS_ABSOLUTE,
        .taxisunit = -1,
        .suffix = "grb",
        .prefix = default_fname_prefix,
      };
#endif

static const struct
{
  char suffix[5];
  int type, defaultDT, defaultGrid;
} suffix2type[] = {
  { "nc", CDI_FILETYPE_NC, CDI_DATATYPE_FLT64, GRID_LONLAT },
  { "grb", CDI_FILETYPE_GRB, CDI_DATATYPE_PACK24, GRID_LONLAT },
  { "grb2", CDI_FILETYPE_GRB2, CDI_DATATYPE_PACK24, GRID_LONLAT },
  { "nc2", CDI_FILETYPE_NC2, CDI_DATATYPE_FLT64, GRID_LONLAT },
  { "nc4", CDI_FILETYPE_NC4, CDI_DATATYPE_FLT64, GRID_LONLAT },
  {
      "ext",
      CDI_FILETYPE_EXT,
      CDI_DATATYPE_FLT64,
      GRID_GENERIC,
  },
  {
      "srv",
      CDI_FILETYPE_SRV,
      CDI_DATATYPE_FLT64,
      GRID_GENERIC,
  },
  { "ieg", CDI_FILETYPE_IEG, CDI_DATATYPE_FLT64, GRID_LONLAT },
};

static void invalidOptionDie(const char *format, ...) __attribute__((noreturn));

static void
invalidOptionDie(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  exit(EXIT_FAILURE);
}

static int
parse_intarg(const char msg[])
{
  char *end;
  long temp = strtol(optarg, &end, 0);
  if ((errno == ERANGE && (temp == LONG_MAX || temp == LONG_MIN)) || (errno != 0 && temp == 0))
    {
      perror(msg);
      exit(EXIT_FAILURE);
    }
  if (temp > INT_MAX || temp < INT_MIN) invalidOptionDie("range error: %ld\n", temp);
  return (int) temp;
}

static unsigned
parse_unsignedarg(const char msg[])
{
  char *end;
  unsigned long temp = strtoul(optarg, &end, 0);
  if ((errno == ERANGE && (temp == ULONG_MAX)) || (errno != 0 && temp == 0))
    {
      perror(msg);
      exit(EXIT_FAILURE);
    }
  if (temp > UINT_MAX) invalidOptionDie("range error: %ld\n", temp);
  return (unsigned) temp;
}

typedef int (*pioRoleFunc)(MPI_Comm commSuper, int IOMode, int nProcsIO);

struct boolOptionParse
{
  bool matched, value, invert;
};

static struct boolOptionParse
parseBooleanLongOption(size_t optionStrSize, const char *optionStr, const char *argStr)
{
  bool invert;
  struct boolOptionParse result;
  if ((invert = !strncmp(argStr, optionStr, optionStrSize - 1)) || !strncmp(argStr, optionStr + 3, optionStrSize - 4))
    {
      result.matched = true;
      result.invert = invert;
      size_t valOfs = optionStrSize - (invert ? 1 : 4);
      if (argStr[valOfs] == '\0' || !strcmp(argStr + valOfs, "=true"))
        result.value = !invert;
      else if (!strcmp(argStr + valOfs, "=false"))
        result.value = invert;
      else
        invalidOptionDie("invalid option argument to -q%s: %s\n", optionStr + (invert ? 0 : 3), argStr + valOfs);
    }
  else
    result = (struct boolOptionParse){ .matched = false, .value = false };
  return result;
}

struct string2int
{
  const char *key;
  int val;
};

static bool
parseLongOptionArgStringToInt(const char *str, size_t optionStrLen, const char *optionStr, size_t mapSize,
                              const struct string2int *map, int argc, char *argv[], int *out)
{
  if (!strncmp(str, optionStr, optionStrLen))
    {
      const char *optargArg;
      if (str[optionStrLen] == '=')
        optargArg = str + optionStrLen + 1;
      else if (str[sizeof(optionStr) - 1] == 0 && optind < argc)
        optargArg = argv[optind++];
      else
        invalidOptionDie("missing argument to -q%s", optionStr);
      for (size_t i = 0; i < mapSize; ++i)
        if (!strcmp(optargArg, map[i].key))
          {
            *out = map[i].val;
            return true;
          }
      invalidOptionDie("unknown argument \"%s\" to -q%s", optargArg, optionStr);
    }
  return false;
}

static void
parse_long_option(struct model_config *restrict setup, int pioConfHandle, pioRoleFunc *pioRoleAssign, const char *str, int argc,
                  char *argv[])
{
#ifndef USE_MPI
  (void) pioConfHandle;
  (void) pioRoleAssign;
#endif
  static const char cacheRedistStr[] = "no-cache-redists", pioRoleSchemeOptionStr[] = "pio-role-scheme",
                    curvilinearGridOptionStr[] = "no-create-curvilinear-grid", uuidCreateOptionStr[] = "no-create-uuid",
                    useDistGridOptionStr[] = "no-use-dist-grid", batchedRmaOptionStr[] = "no-batch-rma",
                    presetDecoOptionStr[] = "no-preset-decomposition", datatypeOptionStr[] = "datatype",
                    prefixOptionStr[] = "prefix", taxistypeOptionStr[] = "taxis-type", taxisunitOptionStr[] = "taxis-unit";
  static const struct string2int datatypeArgMap[] = {
    { "pack", CDI_DATATYPE_PACK },     { "pack1", CDI_DATATYPE_PACK1 },   { "pack2", CDI_DATATYPE_PACK2 },
    { "pack3", CDI_DATATYPE_PACK3 },   { "pack4", CDI_DATATYPE_PACK4 },   { "pack5", CDI_DATATYPE_PACK5 },
    { "pack6", CDI_DATATYPE_PACK6 },   { "pack7", CDI_DATATYPE_PACK7 },   { "pack8", CDI_DATATYPE_PACK8 },
    { "pack9", CDI_DATATYPE_PACK9 },   { "pack10", CDI_DATATYPE_PACK10 }, { "pack11", CDI_DATATYPE_PACK11 },
    { "pack12", CDI_DATATYPE_PACK12 }, { "pack13", CDI_DATATYPE_PACK13 }, { "pack14", CDI_DATATYPE_PACK14 },
    { "pack15", CDI_DATATYPE_PACK15 }, { "pack16", CDI_DATATYPE_PACK16 }, { "pack17", CDI_DATATYPE_PACK17 },
    { "pack18", CDI_DATATYPE_PACK18 }, { "pack19", CDI_DATATYPE_PACK19 }, { "pack20", CDI_DATATYPE_PACK20 },
    { "pack21", CDI_DATATYPE_PACK21 }, { "pack22", CDI_DATATYPE_PACK22 }, { "pack23", CDI_DATATYPE_PACK23 },
    { "pack24", CDI_DATATYPE_PACK24 }, { "pack25", CDI_DATATYPE_PACK25 }, { "pack26", CDI_DATATYPE_PACK26 },
    { "pack27", CDI_DATATYPE_PACK27 }, { "pack28", CDI_DATATYPE_PACK28 }, { "pack29", CDI_DATATYPE_PACK29 },
    { "pack30", CDI_DATATYPE_PACK30 }, { "pack31", CDI_DATATYPE_PACK31 }, { "pack32", CDI_DATATYPE_PACK32 },
    { "cpx32", CDI_DATATYPE_CPX32 },   { "cpx64", CDI_DATATYPE_CPX64 },   { "flt32", CDI_DATATYPE_FLT32 },
    { "flt64", CDI_DATATYPE_FLT64 },   { "int8", CDI_DATATYPE_INT8 },     { "int16", CDI_DATATYPE_INT16 },
    { "int32", CDI_DATATYPE_INT32 },   { "uint8", CDI_DATATYPE_UINT8 },   { "uint16", CDI_DATATYPE_UINT16 },
    { "uint32", CDI_DATATYPE_UINT32 },
  };
  static const struct string2int taxistypeArgMap[] = {
    { "absolute", TAXIS_ABSOLUTE },
    { "relative", TAXIS_RELATIVE },
    { "forecast", TAXIS_FORECAST },
  };
  static const struct string2int taxisunitArgMap[] = {
    { "second", TUNIT_SECOND }, { "minute", TUNIT_MINUTE }, { "quarter", TUNIT_QUARTER }, { "30minutes", TUNIT_30MINUTES },
    { "hour", TUNIT_HOUR },     { "3hours", TUNIT_3HOURS }, { "6hours", TUNIT_6HOURS },   { "12hours", TUNIT_12HOURS },
    { "day", TUNIT_DAY },       { "month", TUNIT_MONTH },   { "year", TUNIT_YEAR },
  };
  enum
  {
    datatypeArgMapSize = sizeof(datatypeArgMap) / sizeof(datatypeArgMap[0]),
    taxistypeArgMapSize = sizeof(taxistypeArgMap) / sizeof(taxistypeArgMap[0]),
    taxisunitArgMapSize = sizeof(taxisunitArgMap) / sizeof(taxisunitArgMap[0]),
  };
  struct boolOptionParse bop;
  if ((bop = parseBooleanLongOption(sizeof(cacheRedistStr), cacheRedistStr, str)).matched)
    {
#ifdef USE_MPI
      cdiPioConfSetRedistCache(pioConfHandle, bop.value);
#else
      invalidOptionDie("CDI-PIO option -q%s unavailable in non-MPI mode\n", cacheRedistStr + (bop.invert ? 0 : 3));
#endif
    }
  else if (parseLongOptionArgStringToInt(str, sizeof(taxistypeOptionStr) - 1, taxistypeOptionStr, taxistypeArgMapSize,
                                         taxistypeArgMap, argc, argv, &setup->taxistype))
    ;
  else if (parseLongOptionArgStringToInt(str, sizeof(taxisunitOptionStr) - 1, taxisunitOptionStr, taxisunitArgMapSize,
                                         taxisunitArgMap, argc, argv, &setup->taxisunit))
    ;
  else if (parseLongOptionArgStringToInt(str, sizeof(datatypeOptionStr) - 1, datatypeOptionStr, datatypeArgMapSize, datatypeArgMap,
                                         argc, argv, &setup->datatype))
    ;
  else if (!strncmp(str, pioRoleSchemeOptionStr, sizeof(pioRoleSchemeOptionStr) - 1))
    {
#ifdef USE_MPI
      static const char pioRoleSchemeLastN[] = "last", pioRoleSchemeFirstN[] = "first", pioRoleSchemeBalanced[] = "balanced";
      if (str[sizeof(pioRoleSchemeOptionStr) - 1] == '=')
        {
          const char *optargArg = str + sizeof(pioRoleSchemeOptionStr);
          if (!strcmp(optargArg, pioRoleSchemeLastN))
            *pioRoleAssign = cdiPioCSRLastN;
          else if (!strcmp(optargArg, pioRoleSchemeFirstN))
            *pioRoleAssign = cdiPioCSRFirstN;
          else if (!strcmp(optargArg, pioRoleSchemeBalanced))
            *pioRoleAssign = cdiPioCSRBalanced;
          else
            invalidOptionDie("unknown scheme argument \"%s\" to -q%s", optargArg, pioRoleSchemeOptionStr);
        }
      else
        invalidOptionDie("long option %s needs argument\n", pioRoleSchemeOptionStr);
#else
      invalidOptionDie("CDI-PIO option -q%s ignored in non-MPI mode\n", pioRoleSchemeOptionStr);
#endif
    }
  else if (!strncmp(str, prefixOptionStr, sizeof(prefixOptionStr) - 1))
    {
      if (str[sizeof(prefixOptionStr) - 1] == '=')
        {
          setup->prefix = str + sizeof(prefixOptionStr);
        }
      else
        invalidOptionDie("long option %s needs argument\n", prefixOptionStr);
    }
  else if ((bop = parseBooleanLongOption(sizeof(useDistGridOptionStr), useDistGridOptionStr, str)).matched)
    {
#if defined USE_MPI && defined HAVE_PPM_DIST_ARRAY_H
      setup->flags = (setup->flags & ~PIO_WRITE_CONFIG_USE_DIST_GRID_FLAG) | (bop.value << PIO_WRITE_CONFIG_USE_DIST_GRID_BIT);
#else
      invalidOptionDie("CDI-PIO option -q%s feature unavailable %s\n", useDistGridOptionStr + (bop.invert ? 0 : 3),
#ifndef USE_MPI
                       "in non-MPI mode"
#else
                       "without PPM feature distributed multi-array"
#endif
      );
#endif
    }
  else if ((bop = parseBooleanLongOption(sizeof(uuidCreateOptionStr), uuidCreateOptionStr, str)).matched)
    {
      setup->flags = (setup->flags & ~PIO_WRITE_CONFIG_CREATE_UUID_FLAG) | (bop.value << PIO_WRITE_CONFIG_CREATE_UUID_BIT);
    }
  else if ((bop = parseBooleanLongOption(sizeof(presetDecoOptionStr), presetDecoOptionStr, str)).matched)
    {
#ifdef USE_MPI
      setup->flags
          = (setup->flags & ~PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_FLAG) | (bop.value << PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_BIT);
#else
      invalidOptionDie("CDI-PIO option -q%s unavailable in non-MPI mode\n", presetDecoOptionStr + (bop.invert ? 0 : 3));
#endif
    }
  else if ((bop = parseBooleanLongOption(sizeof(batchedRmaOptionStr), batchedRmaOptionStr, str)).matched)
    {
#ifdef USE_MPI
      cdiPioConfSetBatchedRMA(pioConfHandle, bop.value);
#else
      invalidOptionDie("CDI-PIO option -q%s unavailable in non-MPI mode\n", batchedRmaOptionStr + (bop.invert ? 0 : 3));
#endif
    }
  else if ((bop = parseBooleanLongOption(sizeof(curvilinearGridOptionStr), curvilinearGridOptionStr, str)).matched)
    {
      setup->flags = (setup->flags & ~PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_FLAG)
                     | (bop.value << PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_BIT);
    }
  else
    invalidOptionDie("unknown long option: %s\n", str);
}

int
main(int argc, char *argv[])
{
  struct model_config setup = default_setup;

  MPI_Comm commModel;
  int pioConfHandle = 0;
  pioRoleFunc pioRoleAssign = 0;
#ifdef USE_MPI
  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int IOMode = PIO_MPI;
  int nProcsIO = 2;

  xmpi(MPI_Init(&argc, &argv));
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi(MPI_Comm_set_errhandler(commGlob, MPI_ERRORS_RETURN));
  xmpi(MPI_Comm_size(commGlob, &sizeGlob));
  xmpi(MPI_Comm_rank(commGlob, &rankGlob));

  pioConfHandle = cdiPioConfCreate();
  pioRoleAssign = cdiPioCSRLastN;
#endif

  /* seed random generator */
  unsigned randSeed;
  {
#ifdef USE_MPI
    if (rankGlob == 0)
#endif
      {
        struct timeval tv;
        int status = gettimeofday(&tv, NULL);
        if (status != 0)
          {
            perror("failed seed generation!");
            exit(1);
          }
        randSeed = (unsigned) (tv.tv_sec ^ tv.tv_usec);
      }
  }

  {
    int opt;
    while ((opt = getopt(argc, argv,
                         "f:m:n:z:t:y:cs:q:"
#ifdef USE_MPI
                         "p:w:"
#endif
                         ))
           != -1)
      switch (opt)
        {
#ifdef USE_MPI
        case 'p':
          IOMode = cdiPioStr2IOMode(optarg);
          if (IOMode < 0)
            {
              fprintf(stderr, "Unsupported PIO mode requested: %s\n", optarg);
              exit(EXIT_FAILURE);
            }
          break;
        case 'w':
          {
            long temp = strtol(optarg, NULL, 0);
            if (temp < 0 || temp > INT_MAX / 2)
              {
                fprintf(stderr, "Unsupported number of I/O servers: %ld\n", temp);
                exit(EXIT_FAILURE);
              }
            nProcsIO = (int) temp;
          }
          break;
#endif
        case 'q': parse_long_option(&setup, pioConfHandle, &pioRoleAssign, optarg, argc, argv); break;
        case 'f':
          {
            int found = 0;
            for (size_t i = 0; i < sizeof(suffix2type) / sizeof(suffix2type[0]); ++i)
              if (!strcmp(optarg, suffix2type[i].suffix))
                {
                  found = 1;
                  setup.filetype = suffix2type[i].type;
                  setup.suffix = suffix2type[i].suffix;
                  setup.datatype = suffix2type[i].defaultDT;
                  break;
                }
            if (!found)
              {
                fprintf(stderr, "Unsupported format requested: %s\n", optarg);
                exit(EXIT_FAILURE);
              }
          }
          break;
        case 'm': setup.nlon = parse_intarg("error parsing number of longitudes"); break;
        case 'n': setup.nlat = parse_intarg("error parsing number of latitudes"); break;
        case 'y':
          setup.nvars = parse_intarg("error parsing number of variables");
          if (setup.nvars < 1)
            {
              fputs("number of variables must be greater than zero!\n", stderr);
              exit(EXIT_FAILURE);
            }
          if (setup.nvars > 127)
            {
              fputs("number of variables must not exceed 127!\n", stderr);
              exit(EXIT_FAILURE);
            }
          break;
        case 'z':
          setup.max_nlev = parse_intarg("error parsing number of levels");
          if (setup.max_nlev < 1)
            {
              fputs("number of levels must be greater than zero!\n", stderr);
              exit(EXIT_FAILURE);
            }
          break;
        case 't': setup.nts = parse_intarg("error parsing number of timesteps"); break;
        case 'c': setup.flags &= ~PIO_WRITE_CONFIG_CHECKSUM_FLAG; break;
        case 's': randSeed = parse_unsignedarg("error parsing random seed"); break;
        default: /* '?' */
          fprintf(stderr,
                  "Usage: %s "
                  "[-m nlon] [-n nlat] [-z nlev] [-t nts] [-y num_vars]"
#ifdef USE_MPI
                  " [-p PIO_MODE] [-w NIOSERVERS] [-c]"
#endif
                  "\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
  }

#ifdef USE_MPI
  MPI_Bcast(&randSeed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif
  cdi_seed_repeatable_random(randSeed);
#ifdef USE_MPI
  if (rankGlob == 0)
#endif
    fprintf(stderr, "random seed=%u\n", randSeed);

#ifdef USE_MPI
  int pioNamespace;
  cdiPioConfSetIOMode(pioConfHandle, IOMode);
  cdiPioConfSetCSRole(pioConfHandle, pioRoleAssign(commGlob, IOMode, nProcsIO));
  cdiPioConfSetPartInflate(pioConfHandle, 1.0f);
  int initNamespace = namespaceGetActive();
  commModel = cdiPioInit(commGlob, pioConfHandle, &pioNamespace);
  if (commModel != MPI_COMM_NULL)
    {
      namespaceSetActive(pioNamespace);
#else
  commModel = -1;
#endif

      modelRun(&setup, commModel);

#ifdef USE_MPI
    }
  cdi_repeatable_finalize();
  pioFinalize();
  namespaceSetActive(initNamespace);
  cdiPioConfDestroy(pioConfHandle);
  xt_finalize();
  MPI_Finalize();
#endif
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
