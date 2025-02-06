#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "cdi.h"

const char *
cdiStringError(int cdiErrno)
{
  // clang-format off
  static const char UnknownError[] = "Unknown Error";
  static const char _ETMOF[]       = "Too many open files";
  static const char _EINVAL[]      = "Invalid argument";
  static const char _EISDIR[]      = "Is a directory";
  static const char _EISEMPTY[]    = "File is empty";
  static const char _EUFTYPE[]     = "Unsupported file type";
  static const char _ELIBNAVAIL[]  = "Unsupported file type (library support not compiled in)";
  static const char _EUFSTRUCT[]   = "Unsupported file structure";
  static const char _EUNC4[]       = "Unsupported NetCDF4 structure";
  static const char _EDIMSIZE[]    = "Invalid dimension size";
  static const char _EQENF[]       = "Query entries not found";
  static const char _EQNAVAIL[]    = "Query not available for file type";
  static const char _ELIMIT[]      = "Internal limits exceeded";

  switch (cdiErrno) {
  case CDI_ESYSTEM:
    {
      const char *cp = strerror(errno);
      if (cp == NULL) break;
      return cp;
    }
  case CDI_ETMOF:      return _ETMOF;
  case CDI_EINVAL:     return _EINVAL;
  case CDI_EISDIR:     return _EISDIR;
  case CDI_EISEMPTY:   return _EISEMPTY;
  case CDI_EUFTYPE:    return _EUFTYPE;
  case CDI_ELIBNAVAIL: return _ELIBNAVAIL;
  case CDI_EUFSTRUCT:  return _EUFSTRUCT;
  case CDI_EUNC4:      return _EUNC4;
  case CDI_EDIMSIZE:   return _EDIMSIZE;
  case CDI_EQENF:      return _EQENF;
  case CDI_EQNAVAIL:   return _EQNAVAIL;
  case CDI_ELIMIT:     return _ELIMIT;
  }
  // clang-format on
  return UnknownError;
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
