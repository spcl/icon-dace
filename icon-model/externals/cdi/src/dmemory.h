#ifndef _DMEMORY_H
#define _DMEMORY_H

// clang-format off

#include <stdio.h>

// if DEBUG_MEMORY is defined setenv MEMORY_DEBUG to debug memory
#define  DEBUG_MEMORY

#ifndef  WITH_FUNCTION_NAME
#define  WITH_FUNCTION_NAME
#endif

#ifdef  __cplusplus
extern "C" {
#endif

extern size_t  memTotal(void);
extern void    memDebug(int debug);
extern void    memExitOnError(void);

extern void   *memRealloc(void *ptr, size_t size, const char *file, const char *functionname, int line);
extern void   *memCalloc(size_t nobjs, size_t size, const char *file, const char *functionname, int line);
extern void   *memMalloc(size_t size, const char *file, const char *functionname, int line);
extern void    memFree(void *ptr, const char *file, const char *functionname, int line);

#ifdef  __cplusplus
}
#endif

#ifdef  DEBUG_MEMORY

#ifdef  WITH_FUNCTION_NAME
#define  Realloc(p, s)  memRealloc((p), (s), __FILE__, __func__, __LINE__)
#define   Calloc(n, s)   memCalloc((n), (s), __FILE__, __func__, __LINE__)
#define   Malloc(s)      memMalloc((s), __FILE__, __func__, __LINE__)
#define     Free(p)        memFree((p), __FILE__, __func__, __LINE__)
#else
#define  Realloc(p, s)  memRealloc((p), (s), __FILE__, (void *) NULL, __LINE__)
#define   Calloc(n, s)   memCalloc((n), (s), __FILE__, (void *) NULL, __LINE__)
#define   Malloc(s)      memMalloc((s), __FILE__, (void *) NULL, __LINE__)
#define     Free(p)        memFree((p), __FILE__, (void *) NULL, __LINE__)
#endif

#else

#include <stdlib.h>

#define  Realloc(p, s)  realloc((p), (s))
#define   Calloc(n, s)   calloc((n), (s))
#define   Malloc(s)      malloc((s))
#define     Free(p)        free((p))

#endif /* DEBUG_MEMORY */

// clang-format on

#endif /* _DMEMORY_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
