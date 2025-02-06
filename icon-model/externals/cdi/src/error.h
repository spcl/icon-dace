#ifndef _ERROR_H_
#define _ERROR_H_

// clang-format off

#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef  WITH_CALLER_NAME
#define  WITH_CALLER_NAME
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int _ExitOnError;  // If set to 1, exit on error (default 1)
extern int _Verbose;      // If set to 1, errors are reported (default 1)
extern int _Debug;        // If set to 1, debuggig (default 0)

void SysError_(const char *caller, const char *fmt, ...);
void    Error_(const char *caller, const char *fmt, ...);
void  Warning_(const char *caller, const char *fmt, ...);
/* delegate used by Warning_ unless mode is PIO */
void cdiWarning(const char *caller, const char *fmt, va_list ap);
void  Message_(const char *caller, const char *fmt, ...);

#ifdef WITH_CALLER_NAME
#  define  SysError(...)  SysError_(__func__, __VA_ARGS__)
#  define    Errorc(...)     Error_(  caller, __VA_ARGS__)
#  define     Error(...)     Error_(__func__, __VA_ARGS__)
#  define   Warning(...)   Warning_(__func__, __VA_ARGS__)
#  define  Messagec(...)   Message_(  caller, __VA_ARGS__)
#  define   Message(...)   Message_(__func__, __VA_ARGS__)
#else
#  define  SysError(...)  SysError_((void *), __VA_ARGS__)
#  define    Errorc(...)     Error_((void *), __VA_ARGS__)
#  define     Error(...)     Error_((void *), __VA_ARGS__)
#  define   Warning(...)   Warning_((void *), __VA_ARGS__)
#  define  Messagec(...)   Message_((void *), __VA_ARGS__)
#  define   Message(...)   Message_((void *), __VA_ARGS__)
#endif

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif

void cdiAbortC(const char *caller, const char *filename,
               const char *functionname, int line,
               const char *errorString, ... )
  __attribute__((noreturn));
#define xabortC(caller, ...)                                    \
  cdiAbortC(caller, __FILE__, __func__, __LINE__, __VA_ARGS__ )
#define xabort(...)                                             \
  cdiAbortC(NULL, __FILE__, __func__, __LINE__, __VA_ARGS__ )
#define cdiAbort(file, func, line, ...)                 \
  cdiAbortC(NULL, (file), (func), (line), __VA_ARGS__)

#define xassert(arg) do {                       \
    if ((arg)) { } else {                       \
      xabort("assertion `" #arg "` failed");}   \
  } while(0)

void
cdiAbortC_serial(const char *caller, const char *filename,
                 const char *functionname, int line,
                 const char *errorString, va_list ap)
  __attribute__((noreturn));


bool cdiObsoleteInfo(const char *oldFunction, const char *newFunction);

#if defined (__cplusplus)
}
#endif

// clang-format on

#endif /* _ERROR_H_ */
