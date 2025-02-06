#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>

#include "cdf.h"
#include "cdi.h"
#include "cdi_int.h"
#include "cdf_int.h"

#ifdef HAVE_LIBNETCDF
const char *
cdfLibraryVersion(void)
{
  return nc_inq_libvers();
}

int CDF_Debug = 0;  // If set to 1, debugging

void
cdfDebug(int debug)
{
  CDF_Debug = debug;

  if (CDF_Debug) Message("debug level %d", debug);
}

static void
cdfComment(int ncid)
{
  static char comment[256] = "Climate Data Interface version ";
  static bool init = false;

  if (!init)
    {
      init = true;
      const char *libvers = cdiLibraryVersion();

      if (!isdigit((int) *libvers))
        strcat(comment, "??");
      else
        strcat(comment, libvers);
      strcat(comment, " (https://mpimet.mpg.de/cdi)");
    }

  cdf_put_att_text(ncid, NC_GLOBAL, "CDI", strlen(comment), comment);
}

static bool
has_uri_scheme(const char *uri)
{
  const char *pos = strstr(uri, "://");
  if (pos)
    {
      int len = pos - uri;
      if (strncmp(uri, "file", len) == 0 || strncmp(uri, "https", len) == 0 || strncmp(uri, "s3", len) == 0) return true;
    }

  return false;
}

static int
cdf_open_read(const char *filename, int *filetype)
{
  int ncid = -1;
  int readmode = NC_NOWRITE;
  int status = cdf_open(filename, readmode, &ncid);
  if (status > 0 && ncid < 0) ncid = CDI_ESYSTEM;
#ifdef HAVE_NETCDF4
  else
    {
      int format = -1;
      status = nc_inq_format(ncid, &format);
      if (status == NC_NOERR && format == NC_FORMAT_NETCDF4_CLASSIC) *filetype = CDI_FILETYPE_NC4C;

#ifdef NC_FORMATX_NCZARR
      int modeNC;
      status = nc_inq_format_extended(ncid, &format, &modeNC);
      if (status == NC_NOERR && format == NC_FORMATX_NCZARR) *filetype = CDI_FILETYPE_NCZARR;
#endif
    }
#endif

  return ncid;
}

static int
cdf_open_write(const char *filename, int *filetype)
{
  int ncid = -1;
  int writemode = NC_CLOBBER;

#ifdef NC_64BIT_OFFSET
  if (*filetype == CDI_FILETYPE_NC2) writemode |= NC_64BIT_OFFSET;
#endif
#ifdef NC_64BIT_DATA
  if (*filetype == CDI_FILETYPE_NC5) writemode |= NC_64BIT_DATA;
#endif
#ifdef HAVE_NETCDF4
  if (*filetype == CDI_FILETYPE_NC4C) writemode |= (NC_NETCDF4 | NC_CLASSIC_MODEL);
  if (*filetype == CDI_FILETYPE_NC4) writemode |= NC_NETCDF4;
  if (*filetype == CDI_FILETYPE_NCZARR) writemode |= NC_NETCDF4;
#endif
  if (*filetype == CDI_FILETYPE_NCZARR)
    {
      if (!has_uri_scheme(filename))
        {
          fprintf(stderr, "URI scheme is missing in NCZarr path!\n");
          return CDI_EINVAL;
        }

      cdf_create(filename, writemode, &ncid);
    }
  else
    {
      if (has_uri_scheme(filename)) fprintf(stderr, "URI scheme defined for non NCZarr Data Model!\n");

      cdf__create(filename, writemode, &ncid);
    }

  return ncid;
}

static int
cdfOpenFile(const char *filename, const char *mode, int *filetype)
{
  int ncid = -1;

  if (filename == NULL)
    {
      ncid = CDI_EINVAL;
    }
  else
    {
      int fmode = tolower(*mode);
      switch (fmode)
        {
        case 'r': ncid = cdf_open_read(filename, filetype); break;
        case 'w':
          ncid = cdf_open_write(filename, filetype);
          if (ncid != CDI_EINVAL)
            {
              if (CDI_Version_Info) cdfComment(ncid);
              cdf_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.6");
            }
          break;
        case 'a': cdf_open(filename, NC_WRITE, &ncid); break;
        default: ncid = CDI_EINVAL;
        }
    }

  return ncid;
}

int
cdfOpen(const char *filename, const char *mode, int filetype)
{
  int fileID = -1;
  bool open_file = true;

  if (CDF_Debug) Message("Open %s with mode %c", filename, *mode);

#ifndef NC_64BIT_OFFSET
  if (filetype == CDI_FILETYPE_NC2) open_file = false;
#endif
#ifndef NC_64BIT_DATA
  if (filetype == CDI_FILETYPE_NC5) open_file = false;
#endif

  if (open_file)
    {
      fileID = cdfOpenFile(filename, mode, &filetype);

      if (CDF_Debug) Message("File %s opened with id %d", filename, fileID);
    }
  else
    {
      fileID = CDI_ELIBNAVAIL;
    }

  return fileID;
}

int
cdf4Open(const char *filename, const char *mode, int *filetype)
{
  if (CDF_Debug) Message("Open %s with mode %c", filename, *mode);

#ifdef HAVE_NETCDF4
  int fileID = cdfOpenFile(filename, mode, filetype);
  if (CDF_Debug) Message("File %s opened with id %d", filename, fileID);
  return fileID;
#else
  return CDI_ELIBNAVAIL;
#endif
}

static void
cdfCloseFile(int fileID)
{
  cdf_close(fileID);
}

void
cdfClose(int fileID)
{
  cdfCloseFile(fileID);
}
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
