#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <sys/stat.h>

#include "cdi.h"
#include "cdi_int.h"
#include "cdf.h"
#include "cdf_int.h"
#include "namespace.h"

#ifdef HAVE_LIBNETCDF

void
cdf_create(const char *path, int cmode, int *ncidp)
{
  int status = nc_create(path, cmode, ncidp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  mode=%d  file=%s", *ncidp, cmode, path);

  if (status != NC_NOERR) Error("%s: %s", path, nc_strerror(status));

  int oldfill;
  status = nc_set_fill(*ncidp, NC_NOFILL, &oldfill);

  if (status != NC_NOERR) Error("%s: %s", path, nc_strerror(status));
}

void
cdf__create(const char *path, int cmode, int *ncidp)
{
  int status = -1;
  size_t chunksizehint = 0;

  size_t initialsz = 0;

#if defined(__SX__) || defined(ES)
  chunksizehint = 16777216;  // 16 MB
#endif

  if (CDI_Netcdf_Chunksizehint != CDI_UNDEFID) chunksizehint = (size_t) CDI_Netcdf_Chunksizehint;

  cdi_nc__create_funcp my_nc__create = (cdi_nc__create_funcp) namespaceSwitchGet(NSSWITCH_NC__CREATE).func;
  status = my_nc__create(path, cmode, initialsz, &chunksizehint, ncidp);

  if (status != NC_NOERR)
    {
      if (CDF_Debug) Message("ncid=%d  mode=%d  chunksizehint=%zu  file=%s", *ncidp, cmode, chunksizehint, path);
      Error("%s: %s", path, nc_strerror(status));
    }

  int oldfill;
  status = nc_set_fill(*ncidp, NC_NOFILL, &oldfill);

  if (status != NC_NOERR) Error("%s: %s", path, nc_strerror(status));
}

int
cdf_open(const char *path, int omode, int *ncidp)
{
  int status = 0;

  if (strstr(path, ":/"))  // ESDM and DAP
    {
      status = nc_open(path, omode, ncidp);
    }
  else
    {
      struct stat filestat;
      if (stat(path, &filestat) != 0) SysError(path);

      size_t chunksizehint = 0;
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
      chunksizehint = (size_t) filestat.st_blksize * 4;
      if (chunksizehint > (size_t) filestat.st_size) chunksizehint = (size_t) filestat.st_size;
#endif
      // if (chunksizehint < ChunkSizeMin) chunksizehint = ChunkSizeMin;
      if (CDI_Netcdf_Chunksizehint != CDI_UNDEFID) chunksizehint = (size_t) CDI_Netcdf_Chunksizehint;

      // FIXME: parallel part missing
      status = nc__open(path, omode, &chunksizehint, ncidp);

      if (CDF_Debug) Message("chunksizehint %zu", chunksizehint);
    }

  if (CDF_Debug) Message("ncid=%d  mode=%d  file=%s", *ncidp, omode, path);

  if (CDF_Debug && status != NC_NOERR) Message("%s", nc_strerror(status));

  return status;
}

void
cdf_close(int ncid)
{
  int status = nc_close(ncid);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_redef(int ncid)
{
  int status = nc_redef(ncid);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

int
cdi_nc_enddef_serial(int ncid, int streamID)
{
  (void) streamID;
  return nc_enddef(ncid);
}

int
cdi_nc__enddef_serial(int ncid, int streamID, size_t hdr_pad, size_t v_align, size_t v_minfree, size_t r_align)
{
  (void) streamID;
  return nc__enddef(ncid, hdr_pad, v_align, v_minfree, r_align);
}

void
cdf_enddef(int ncid, int streamID)
{
  cdi_nc_enddef_funcp my_nc_enddef = (cdi_nc_enddef_funcp) namespaceSwitchGet(NSSWITCH_NC_ENDDEF).func;
  int status = my_nc_enddef(ncid, streamID);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf__enddef(int ncid, int streamID, const size_t hdr_pad)
{
  const size_t v_align = 4UL;    // [B] Alignment of beginning of data section for fixed variables
  const size_t v_minfree = 0UL;  // [B] Pad at end of data section for fixed size variables
  const size_t r_align = 4UL;    // [B] Alignment of beginning of data section for record variables

  // nc_enddef(ncid) is equivalent to nc__enddef(ncid, 0, 4, 0, 4)
  cdi_nc__enddef_funcp my_nc__enddef = (cdi_nc__enddef_funcp) namespaceSwitchGet(NSSWITCH_NC_ENDDEF).func;
  int status = my_nc__enddef(ncid, streamID, hdr_pad, v_align, v_minfree, r_align);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_sync(int ncid)
{
  int status = nc_sync(ncid);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq(int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp)
{
  int status = nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);

  if (CDF_Debug || status != NC_NOERR)
    Message("ncid=%d  ndims=%d  nvars=%d  ngatts=%d  unlimid=%d", ncid, *ndimsp, *nvarsp, *ngattsp, *unlimdimidp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_def_dim(int ncid, const char *name, size_t len, int *dimidp)
{
  int status = nc_def_dim(ncid, name, len, dimidp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  name=%s  len=%d", ncid, name, len);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_dimid(int ncid, const char *name, int *dimidp)
{
  int status = nc_inq_dimid(ncid, name, dimidp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  name=%s  dimid=%d", ncid, name, *dimidp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_dim(int ncid, int dimid, char *name, size_t *lengthp)
{
  int status = nc_inq_dim(ncid, dimid, name, lengthp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  dimid=%d  length=%d  name=%s", ncid, dimid, *lengthp, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_dimname(int ncid, int dimid, char *name)
{
  int status = nc_inq_dimname(ncid, dimid, name);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  dimid=%d  name=%s", ncid, dimid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_dimlen(int ncid, int dimid, size_t *lengthp)
{
  int status = nc_inq_dimlen(ncid, dimid, lengthp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  dimid=%d  length=%d", ncid, dimid, *lengthp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_def_var(int ncid, const char *name, nc_type xtype, int ndims, const int dimids[], int *varidp)
{
  cdi_cdf_def_var_funcp my_cdf_def_var = (cdi_cdf_def_var_funcp) namespaceSwitchGet(NSSWITCH_CDF_DEF_VAR).func;
  my_cdf_def_var(ncid, name, xtype, ndims, dimids, varidp);
}

void
cdf_def_var_serial(int ncid, const char *name, nc_type xtype, int ndims, const int dimids[], int *varidp)
{
  int status = nc_def_var(ncid, name, xtype, ndims, dimids, varidp);
  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  name=%s  xtype=%d  ndims=%d  varid=%d", ncid, name, xtype, ndims, *varidp);
  if (status == NC_NOERR)
    {
      int fileFormat;
      status = nc_inq_format(ncid, &fileFormat);
      if (status == NC_NOERR && (fileFormat == NC_FORMAT_NETCDF4 || fileFormat == NC_FORMAT_NETCDF4_CLASSIC))
        status = nc_def_var_fill(ncid, *varidp, 1, NULL);
    }
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_varid(int ncid, const char *name, int *varidp)
{
  int status = nc_inq_varid(ncid, name, varidp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  name=%s  varid=%d", ncid, name, *varidp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_nvars(int ncid, int *nvarsp)
{
  int status = nc_inq_nvars(ncid, nvarsp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  nvars=%d", ncid, *nvarsp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_var(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int dimids[], int *nattsp)
{
  int status = nc_inq_var(ncid, varid, name, xtypep, ndimsp, dimids, nattsp);

  if (CDF_Debug || status != NC_NOERR)
    Message("ncid=%d  varid=%d  ndims=%d  xtype=%d  natts=%d  name=%s", ncid, varid, *ndimsp, *xtypep, *nattsp, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_varname(int ncid, int varid, char *name)
{
  int status = nc_inq_varname(ncid, varid, name);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  name=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_vartype(int ncid, int varid, nc_type *xtypep)
{
  int status = nc_inq_vartype(ncid, varid, xtypep);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  xtype=%s", ncid, varid, *xtypep);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_varndims(int ncid, int varid, int *ndimsp)
{
  int status = nc_inq_varndims(ncid, varid, ndimsp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_vardimid(int ncid, int varid, int dimids[])
{
  int status = nc_inq_vardimid(ncid, varid, dimids);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_varnatts(int ncid, int varid, int *nattsp)
{
  int status = nc_inq_varnatts(ncid, varid, nattsp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  nattsp=%d", ncid, varid, *nattsp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_text(int ncid, int varid, const char *tp)
{
  int status = nc_put_var_text(ncid, varid, tp);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %s", ncid, varid, tp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_short(int ncid, int varid, const short *sp)
{
  int status = nc_put_var_short(ncid, varid, sp);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %hd", ncid, varid, *sp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_int(int ncid, int varid, const int *ip)
{
  int status = nc_put_var_int(ncid, varid, ip);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %d", ncid, varid, *ip);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_long(int ncid, int varid, const long *lp)
{
  int status = nc_put_var_long(ncid, varid, lp);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %ld", ncid, varid, *lp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_float(int ncid, int varid, const float *fp)
{
  int status = nc_put_var_float(ncid, varid, fp);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %f", ncid, varid, *fp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

static const char *
cdf_var_type(nc_type xtype)
{
  const char *ctype = "unknown";

  // clang-format off
  if      (xtype == NC_BYTE  )  ctype = "NC_BYTE";
  else if (xtype == NC_CHAR  )  ctype = "NC_CHAR";
  else if (xtype == NC_SHORT )  ctype = "NC_SHORT";
  else if (xtype == NC_INT   )  ctype = "NC_INT";
  else if (xtype == NC_FLOAT )  ctype = "NC_FLOAT";
  else if (xtype == NC_DOUBLE)  ctype = "NC_DOUBLE";
#ifdef  HAVE_NETCDF4
  else if (xtype == NC_UBYTE )  ctype = "NC_UBYTE";
  else if (xtype == NC_LONG  )  ctype = "NC_LONG";
  else if (xtype == NC_USHORT)  ctype = "NC_USHORT";
  else if (xtype == NC_UINT  )  ctype = "NC_UINT";
  else if (xtype == NC_INT64 )  ctype = "NC_INT64";
  else if (xtype == NC_UINT64)  ctype = "NC_UINT64";
#endif
  // clang-format on

  return ctype;
}

static void
minmaxval(size_t nvals, const double *array, double *minval, double *maxval)
{
  double minv = array[0];
  double maxv = array[0];
  for (size_t i = 0; i < nvals; ++i)
    {
      minv = (array[i] < minv) ? array[i] : minv;
      maxv = (array[i] > maxv) ? array[i] : maxv;
    }

  *minval = minv;
  *maxval = maxv;
}

static void
minmaxvalf(size_t nvals, const float *array, double *minval, double *maxval)
{
  float minv = array[0];
  float maxv = array[0];
  for (size_t i = 0; i < nvals; ++i)
    {
      minv = (array[i] < minv) ? array[i] : minv;
      maxv = (array[i] > maxv) ? array[i] : maxv;
    }

  *minval = minv;
  *maxval = maxv;
}

void
cdf_put_vara_double(int ncid, int varid, const size_t start[], const size_t count[], const double *dp)
{
  int status = nc_put_vara_double(ncid, varid, start, count, dp);

  if (CDF_Debug || status != NC_NOERR)
    {
      char name[256];
      nc_inq_varname(ncid, varid, name);
      nc_type xtype;
      nc_inq_vartype(ncid, varid, &xtype);
      int ndims;
      nc_inq_varndims(ncid, varid, &ndims);
      double minval = 0.0, maxval = 0.0;
      size_t nvals = 1;
      for (int i = 0; i < ndims; ++i) nvals *= count[i];
      minmaxval(nvals, dp, &minval, &maxval);
      Message("name=%s  type=%s  minval=%f  maxval=%f", name, cdf_var_type(xtype), minval, maxval);
    }

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_vara_float(int ncid, int varid, const size_t start[], const size_t count[], const float *fp)
{
  int status = nc_put_vara_float(ncid, varid, start, count, fp);

  if (CDF_Debug || status != NC_NOERR)
    {
      char name[256];
      nc_inq_varname(ncid, varid, name);
      nc_type xtype;
      nc_inq_vartype(ncid, varid, &xtype);
      int ndims;
      nc_inq_varndims(ncid, varid, &ndims);
      double minval = 0.0, maxval = 0.0;
      size_t nvals = 1;
      for (int i = 0; i < ndims; ++i) nvals *= count[i];
      minmaxvalf(nvals, fp, &minval, &maxval);
      Message("name=%s  type=%s  minval=%f  maxval=%f", name, cdf_var_type(xtype), minval, maxval);
    }

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_vara(int ncid, int varid, const size_t start[], const size_t count[], const void *cp)
{
  int status = nc_put_vara(ncid, varid, start, count, cp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara(int ncid, int varid, const size_t start[], const size_t count[], void *cp)
{
  int status = nc_get_vara(ncid, varid, start, count, cp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara_int(int ncid, int varid, const size_t start[], const size_t count[], int *dp)
{
  int status = nc_get_vara_int(ncid, varid, start, count, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara_double(int ncid, int varid, const size_t start[], const size_t count[], double *dp)
{
  int status = nc_get_vara_double(ncid, varid, start, count, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  start[0]=%zu  count[0]=%zu", ncid, varid, start[0], count[0]);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara_float(int ncid, int varid, const size_t start[], const size_t count[], float *fp)
{
  int status = nc_get_vara_float(ncid, varid, start, count, fp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  start[0]=%zu  count[0]=%zu", ncid, varid, start[0], count[0]);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara_text(int ncid, int varid, const size_t start[], const size_t count[], char *tp)
{
  int status = nc_get_vara_text(ncid, varid, start, count, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_vara_uchar(int ncid, int varid, const size_t start[], const size_t count[], unsigned char *tp)
{
  int status = nc_get_vara_uchar(ncid, varid, start, count, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var_double(int ncid, int varid, const double *dp)
{
  int status = nc_put_var_double(ncid, varid, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  val0=%f", ncid, varid, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var1_text(int ncid, int varid, const size_t index[], char *tp)
{
  int status = nc_get_var1_text(ncid, varid, index, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var1_double(int ncid, int varid, const size_t index[], double *dp)
{
  int status = nc_get_var1_double(ncid, varid, index, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_var1_double(int ncid, int varid, const size_t index[], const double *dp)
{
  int status = nc_put_var1_double(ncid, varid, index, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  val=%f", ncid, varid, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_text(int ncid, int varid, char *tp)
{
  int status = nc_get_var_text(ncid, varid, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_short(int ncid, int varid, short *sp)
{
  int status = nc_get_var_short(ncid, varid, sp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_int(int ncid, int varid, int *ip)
{
  int status = nc_get_var_int(ncid, varid, ip);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_long(int ncid, int varid, long *lp)
{
  int status = nc_get_var_long(ncid, varid, lp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_float(int ncid, int varid, float *fp)
{
  int status = nc_get_var_float(ncid, varid, fp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d", ncid, varid);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_var_double(int ncid, int varid, double *dp)
{
  int status = nc_get_var_double(ncid, varid, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  val[0]=%f", ncid, varid, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out, int varid_out)
{
  int status = nc_copy_att(ncid_in, varid_in, name, ncid_out, varid_out);

  if (CDF_Debug || status != NC_NOERR) Message("%d %d %s %d %d", ncid_in, varid_out, name, ncid_out, varid_out);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_att_text(int ncid, int varid, const char *name, size_t len, const char *tp)
{
  int status = nc_put_att_text(ncid, varid, name, len, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  text=%.*s", ncid, varid, name, (int) len, tp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_att_int(int ncid, int varid, const char *name, nc_type xtype, size_t len, const int *ip)
{
  int status = nc_put_att_int(ncid, varid, name, xtype, len, ip);

  if (status == NC_ERANGE) status = NC_NOERR;

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%d", ncid, varid, name, *ip);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_att_float(int ncid, int varid, const char *name, nc_type xtype, size_t len, const float *dp)
{
  int status = nc_put_att_float(ncid, varid, name, xtype, len, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%g", ncid, varid, name, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_put_att_double(int ncid, int varid, const char *name, nc_type xtype, size_t len, const double *dp)
{
  int status = nc_put_att_double(ncid, varid, name, xtype, len, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%g", ncid, varid, name, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_att_text(int ncid, int varid, const char *name, char *tp)
{
  int status = nc_get_att_text(ncid, varid, name, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  name=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_att_string(int ncid, int varid, const char *name, char **tp)
{
#ifdef HAVE_NETCDF4
  int status = nc_get_att_string(ncid, varid, name, tp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  name=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
#endif
}

void
cdf_get_att_int(int ncid, int varid, const char *name, int *ip)
{
  int status = nc_get_att_int(ncid, varid, name, ip);

  if (status == NC_ERANGE) status = NC_NOERR;

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%d", ncid, varid, name, *ip);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_get_att_longlong(int ncid, int varid, const char *name, long long *llp)
{
#ifdef HAVE_NETCDF4
  int status = nc_get_att_longlong(ncid, varid, name, llp);

  if (status == NC_ERANGE) status = NC_NOERR;

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%lld", ncid, varid, name, *llp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
#endif
}

void
cdf_get_att_double(int ncid, int varid, const char *name, double *dp)
{
  int status = nc_get_att_double(ncid, varid, name, dp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  val=%.9g", ncid, varid, name, *dp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_att(int ncid, int varid, const char *name, nc_type *xtypep, size_t *lenp)
{
  int status = nc_inq_att(ncid, varid, name, xtypep, lenp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_atttype(int ncid, int varid, const char *name, nc_type *xtypep)
{
  int status = nc_inq_atttype(ncid, varid, name, xtypep);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_attlen(int ncid, int varid, const char *name, size_t *lenp)
{
  int status = nc_inq_attlen(ncid, varid, name, lenp);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s  len=%d", ncid, varid, name, *lenp);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_attname(int ncid, int varid, int attnum, char *name)
{
  int status = nc_inq_attname(ncid, varid, attnum, name);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  attnum=%d  att=%s", ncid, varid, attnum, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

void
cdf_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
  int status = nc_inq_attid(ncid, varid, name, attnump);

  if (CDF_Debug || status != NC_NOERR) Message("ncid=%d  varid=%d  att=%s", ncid, varid, name);

  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}

#ifdef HAVE_NETCDF4
void
cdf_def_var_chunking(int ncid, int varid, int storage, const size_t *chunksizesp)
{
  int status = nc_def_var_chunking(ncid, varid, storage, chunksizesp);
  if (status != NC_NOERR) Error("%s", nc_strerror(status));
}
#endif

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
