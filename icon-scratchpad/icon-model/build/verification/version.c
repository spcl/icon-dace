#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#if defined(__cplusplus)
extern "C"
{
#  define LANG_PREFIX cxx_
#else
#  define LANG_PREFIX c_
#endif

#define GLUE_helper(x, y) x##y
#define GLUE(x, y) GLUE_helper(x, y)

#define CHARS1(x) ('0' + ((x) % 10))
#define CHARS2(x) ('0' + (((x) / 10) % 10)), CHARS1(x)
#define CHARS3(x) ('0' + (((x) / 100) % 10)), CHARS2(x)
#define CHARS4(x) ('0' + (((x) / 1000) % 10)), CHARS3(x)
#define CHARS5(x) ('0' + (((x) / 10000) % 10)), CHARS4(x)

#define DEFINE_STRING_GETTER_helper(x)                             \
  void pvcs_get_##x(const char **const val, size_t *const val_len) \
  {                                                                \
    *val = pvcs_##x;                                               \
    *val_len = (*val == NULL) ? 0 : strlen(*val);                  \
  }
#define DEFINE_STRING_GETTER(x) DEFINE_STRING_GETTER_helper(x)

#if !defined(COMPILER_VERSION_ONLY)

struct pvcs_info
{
  const char *const name;
  const char *const revision;
  const char *const remote_url;
  const char *const local_branch;
};

static const struct pvcs_info pvcs_infos[] = {
  // icon (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model)
  { .name = "icon",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // ecrad (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad)
  { .name = "ecrad",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // cdi (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/cdi)
  { .name = "cdi",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // mtime (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/mtime)
  { .name = "mtime",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // fortran-support (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/fortran-support)
  { .name = "fortran-support",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // math-support (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/math-support)
  { .name = "math-support",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
  // math-interpolation (/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/math-interpolation)
  { .name = "math-interpolation",
    .revision = NULL,
    .remote_url = NULL,
    .local_branch = NULL, },
};

static const struct pvcs_info *
pvcs_find_info(const char *key, const size_t key_len)
{
  static const size_t count = sizeof(pvcs_infos) / sizeof(struct pvcs_info);
  for (size_t ii = 0; ii < count; ++ii)
    {
      const struct pvcs_info *info = &pvcs_infos[ii];
      if (strlen(info->name) == key_len
          && strncmp(info->name, key, key_len) == 0)
        {
          return info;
        }
    }
  return NULL;
}

#  define DEFINE_INFO_GETTER(x)                                          \
    int pvcs_get_##x(const char *const key, const size_t key_len,        \
                     const char **const val, size_t *const val_len)      \
    {                                                                    \
      const struct pvcs_info *const info = pvcs_find_info(key, key_len); \
      *val = (info == NULL) ? NULL : info->x;                            \
      *val_len = (*val == NULL) ? 0 : strlen(info->x);                   \
      return (info == NULL) ? 1 : 0;                                     \
    }

DEFINE_INFO_GETTER(revision)
DEFINE_INFO_GETTER(remote_url)
DEFINE_INFO_GETTER(local_branch)

static const char *const pvcs_icon_version =
#  if defined(PACKAGE_VERSION)
  PACKAGE_VERSION
#  else
  NULL
#  endif
  ;

DEFINE_STRING_GETTER(icon_version)

#endif

// Intel
#if defined(__INTEL_LLVM_COMPILER)
#  define COMPILER_NAME "Intel"
#  if __INTEL_LLVM_COMPILER < 1000000L
#    define COMPILER_VERSION_MAJOR __INTEL_LLVM_COMPILER / 100
#    define COMPILER_VERSION_MINOR __INTEL_LLVM_COMPILER / 10 % 10
#    undef COMPILER_VERSION_PATCH
#  else
#    define COMPILER_VERSION_MAJOR __INTEL_LLVM_COMPILER / 10000
#    define COMPILER_VERSION_MINOR __INTEL_LLVM_COMPILER / 100 % 100
#    define COMPILER_VERSION_PATCH __INTEL_LLVM_COMPILER % 100
#  endif
// Intel Classic
#elif defined(__INTEL_COMPILER)
#  define COMPILER_NAME "Intel Classic"
#  if __INTEL_COMPILER < 2021
#    define COMPILER_VERSION_MAJOR __INTEL_COMPILER / 100
#    define COMPILER_VERSION_MINOR __INTEL_COMPILER / 10 % 10
#    if __INTEL_COMPILER_BUILD_DATE == 20181018 \
      || __INTEL_COMPILER_BUILD_DATE == 20200306
#      define COMPILER_VERSION_PATCH 1
#    elif defined(__INTEL_COMPILER_UPDATE)
#      define COMPILER_VERSION_PATCH __INTEL_COMPILER_UPDATE
#    else
#      undef COMPILER_VERSION_PATCH
#    endif
#  else
#    define COMPILER_VERSION_MAJOR __INTEL_COMPILER
#    define COMPILER_VERSION_MINOR __INTEL_COMPILER_UPDATE
#    if __INTEL_COMPILER_BUILD_DATE == 20201208
#      define COMPILER_VERSION_PATCH 2
#    else
#      define COMPILER_VERSION_PATCH 0
#    endif
#  endif
// Cray
#elif defined(__cray__)
#  define COMPILER_NAME "Cray"
#  define COMPILER_VERSION_MAJOR __cray_major__
#  define COMPILER_VERSION_MINOR __cray_minor__
#  define COMPILER_VERSION_PATCH __cray_patchlevel__
// Cray Classic
#elif defined(_CRAYC)
#  define COMPILER_NAME "Cray Classic"
#  define COMPILER_VERSION_MAJOR _RELEASE_MAJOR
#  define COMPILER_VERSION_MINOR _RELEASE_MINOR
#  define COMPILER_VERSION_PATCH _RELEASE_PATCHLEVEL
// NEC
#elif defined(__NEC__)
#  define COMPILER_NAME "NEC"
#  define COMPILER_VERSION_MAJOR __NEC_VERSION__ / 10000
#  define COMPILER_VERSION_MINOR __NEC_VERSION__ / 100 % 100
#  define COMPILER_VERSION_PATCH __NEC_VERSION__ % 100
// NVHPC
#elif defined(__NVCOMPILER)
#  define COMPILER_NAME "NVHPC"
#  define COMPILER_VERSION_MAJOR __NVCOMPILER_MAJOR__
#  define COMPILER_VERSION_MINOR __NVCOMPILER_MINOR__
#  define COMPILER_VERSION_PATCH __NVCOMPILER_PATCHLEVEL__
// PGI
#elif defined(__PGI)
#  define COMPILER_NAME "PGI"
#  define COMPILER_VERSION_MAJOR __PGIC__
#  define COMPILER_VERSION_MINOR __PGIC_MINOR__
#  define COMPILER_VERSION_PATCH __PGIC_PATCHLEVEL__
// Apple Clang
#elif defined(__apple_build_version__)
#  define COMPILER_NAME "Apple Clang"
#  define COMPILER_VERSION_MAJOR __clang_major__
#  define COMPILER_VERSION_MINOR __clang_minor__
#  define COMPILER_VERSION_PATCH __clang_patchlevel__
// Clang
#elif defined(__clang__)
#  define COMPILER_NAME "Clang"
#  define COMPILER_VERSION_MAJOR __clang_major__
#  define COMPILER_VERSION_MINOR __clang_minor__
#  define COMPILER_VERSION_PATCH __clang_patchlevel__
#else
#  undef COMPILER_NAME
#  undef COMPILER_VERSION_MAJOR
#  undef COMPILER_VERSION_MINOR
#  undef COMPILER_VERSION_PATCH
#endif

// GCC: either the primary or the secondary compiler
#if defined(__GNUC__)
#  if defined(COMPILER_NAME)
#    define COMPILER_SECONDARY_NAME "GNU"
#    define COMPILER_SECONDARY_VERSION_MAJOR __GNUC__
#    define COMPILER_SECONDARY_VERSION_MINOR __GNUC_MINOR__
#    define COMPILER_SECONDARY_VERSION_PATCH __GNUC_PATCHLEVEL__
#  else
#    define COMPILER_NAME "GNU"
#    define COMPILER_VERSION_MAJOR __GNUC__
#    define COMPILER_VERSION_MINOR __GNUC_MINOR__
#    define COMPILER_VERSION_PATCH __GNUC_PATCHLEVEL__
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_name)) =
#if defined(COMPILER_NAME)
  COMPILER_NAME
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_name))

#undef PATCH_CHARS
#undef MINOR_CHARS
#undef MAJOR_CHARS
#if defined(COMPILER_VERSION_MAJOR)
#  if COMPILER_VERSION_MAJOR < 10
#    define MAJOR_CHARS CHARS1(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 100
#    define MAJOR_CHARS CHARS2(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 1000
#    define MAJOR_CHARS CHARS3(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 10000
#    define MAJOR_CHARS CHARS4(COMPILER_VERSION_MAJOR)
#  else
#    define MAJOR_CHARS CHARS5(COMPILER_VERSION_MAJOR)
#  endif
#  if defined(COMPILER_VERSION_MINOR)
#    if COMPILER_VERSION_MINOR < 10
#      define MINOR_CHARS CHARS1(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 100
#      define MINOR_CHARS CHARS2(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 1000
#      define MINOR_CHARS CHARS3(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 10000
#      define MINOR_CHARS CHARS4(COMPILER_VERSION_MINOR)
#    else
#      define MINOR_CHARS CHARS5(COMPILER_VERSION_MINOR)
#    endif
#    if defined(COMPILER_VERSION_PATCH)
#      if COMPILER_VERSION_PATCH < 10
#        define PATCH_CHARS CHARS1(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 100
#        define PATCH_CHARS CHARS2(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 1000
#        define PATCH_CHARS CHARS3(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 10000
#        define PATCH_CHARS CHARS4(COMPILER_VERSION_PATCH)
#      else
#        define PATCH_CHARS CHARS5(COMPILER_VERSION_PATCH)
#      endif
#    endif
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_version)) =
#if defined(MAJOR_CHARS)
  (const char[])
  {
    MAJOR_CHARS,
#  if defined(MINOR_CHARS)
    '.', MINOR_CHARS,
#    if defined(PATCH_CHARS)
    '.', PATCH_CHARS,
#    endif
#  endif
    '\0',
  }
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_version))

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_secondary_name)) =
#if defined(COMPILER_SECONDARY_NAME)
  COMPILER_SECONDARY_NAME
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_secondary_name))

#undef PATCH_CHARS
#undef MINOR_CHARS
#undef MAJOR_CHARS
#if defined(COMPILER_SECONDARY_VERSION_MAJOR)
#  if COMPILER_SECONDARY_VERSION_MAJOR < 10
#    define MAJOR_CHARS CHARS1(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 100
#    define MAJOR_CHARS CHARS2(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 1000
#    define MAJOR_CHARS CHARS3(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 10000
#    define MAJOR_CHARS CHARS4(COMPILER_SECONDARY_VERSION_MAJOR)
#  else
#    define MAJOR_CHARS CHARS5(COMPILER_SECONDARY_VERSION_MAJOR)
#  endif
#  if defined(COMPILER_SECONDARY_VERSION_MINOR)
#    if COMPILER_SECONDARY_VERSION_MINOR < 10
#      define MINOR_CHARS CHARS1(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 100
#      define MINOR_CHARS CHARS2(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 1000
#      define MINOR_CHARS CHARS3(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 10000
#      define MINOR_CHARS CHARS4(COMPILER_SECONDARY_VERSION_MINOR)
#    else
#      define MINOR_CHARS CHARS5(COMPILER_SECONDARY_VERSION_MINOR)
#    endif
#    if defined(COMPILER_SECONDARY_VERSION_PATCH)
#      if COMPILER_SECONDARY_VERSION_PATCH < 10
#        define PATCH_CHARS CHARS1(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 100
#        define PATCH_CHARS CHARS2(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 1000
#        define PATCH_CHARS CHARS3(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 10000
#        define PATCH_CHARS CHARS4(COMPILER_SECONDARY_VERSION_PATCH)
#      else
#        define PATCH_CHARS CHARS5(COMPILER_SECONDARY_VERSION_PATCH)
#      endif
#    endif
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_secondary_version)) =
#if defined(MAJOR_CHARS)
  (const char[])
  {
    MAJOR_CHARS,
#  if defined(MINOR_CHARS)
    '.', MINOR_CHARS,
#    if defined(PATCH_CHARS)
    '.', PATCH_CHARS,
#    endif
#  endif
    '\0',
  }
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_secondary_version))

#if defined(__cplusplus)
}
#endif
