/*
 * @file ppm_posix_c.c
 * @brief Fortran wrappers for POSIX C functions, C part
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords:
 * @author Thomas Jahns <jahns@dkrz.de>
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "core/ppm_visibility.h"
#define FCALLSC_QUALIFIER PPM_DSO_INTERNAL
#include "cfortran.h"
#include "core/minmax.h"

static void
PPM_mkdir(const char *path, int mode, int *ierr)
{
  *ierr = mkdir(path, (mode_t)mode) == 0 ? 0 : errno;
}

FCALLSCSUB3(PPM_mkdir,PPM_MKDIR_F,ppm_mkdir_f,STRING,INT,PINT)

static void
PPM_rmdir(const char *path, int *ierr)
{
  *ierr = rmdir(path) == 0 ? 0 : errno;
}

FCALLSCSUB2(PPM_rmdir,PPM_RMDIR_F,ppm_rmdir_f,STRING,PINT)

struct PPM_stat_f
{
  int32_t ppm_st_dev, ppm_st_mode;
  int64_t ppm_st_ino;
  int32_t ppm_st_nlink, ppm_st_uid, ppm_st_gid, ppm_st_rdev;
  int64_t ppm_st_size;
  int32_t ppm_st_blksize, fill;
  int64_t ppm_st_blocks,
    ppm_st_atime,
    ppm_st_mtime,
    ppm_st_ctime;
};

static void
PPM_stat(const char *path, struct PPM_stat_f *buf, int *ierr)
{
  struct stat buf_temp;
  *ierr = stat(path, &buf_temp) == 0 ? 0 : errno;
  buf->ppm_st_dev = (int32_t)buf_temp.st_dev;
  buf->ppm_st_ino = (int64_t)buf_temp.st_ino;
  buf->ppm_st_mode = (int32_t)buf_temp.st_mode;
  buf->ppm_st_nlink = (int32_t)buf_temp.st_nlink;
  buf->ppm_st_uid = (int32_t)buf_temp.st_uid;
  buf->ppm_st_gid = (int32_t)buf_temp.st_gid;
  buf->ppm_st_rdev = (int32_t)buf_temp.st_rdev;
  buf->ppm_st_size = (int64_t)buf_temp.st_size;
  buf->ppm_st_blksize = (int32_t)buf_temp.st_blksize;
  buf->ppm_st_blocks = (int64_t)buf_temp.st_blocks;
  buf->ppm_st_atime = (int64_t)buf_temp.st_atime;
  buf->ppm_st_mtime = (int64_t)buf_temp.st_mtime;
  buf->ppm_st_ctime = (int64_t)buf_temp.st_ctime;
}

FCALLSCSUB3(PPM_stat,PPM_STAT_F,ppm_stat_f,STRING,PVOID,PINT)

static int
PPM_is_dir(struct PPM_stat_f *stats)
{
  return S_ISDIR((mode_t)stats->ppm_st_mode);
}

FCALLSCFUN1(LOGICAL,PPM_is_dir,PPM_IS_DIR_F,ppm_is_dir_f,PVOID)

static void
PPM_strerror(char *buf, int buf_len, int ierr, int *result_len)
{
  const char *s = ierr != 0 ? strerror(ierr) : "";
  size_t n = MIN3(strlen(s), (size_t)SIZE_MAX, (unsigned)INT_MAX);
  *result_len = (int)n;
  buf_len = buf_len >= 0 ? buf_len : 0;
  n = n > (size_t)buf_len ? (size_t)buf_len : n;
  memcpy(buf, s, n);
  buf[n] = '\0';
}

FCALLSCSUB4(PPM_strerror,PPM_STRERROR_F,ppm_strerror_f,PSTRING,INT,INT,PINT)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
