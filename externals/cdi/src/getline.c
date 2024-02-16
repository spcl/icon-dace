/*
 * getline.c --- replacement for GNU getline if not available in libc
 *
 * Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
 *
 * Version: 0.0.1
 * Author: Thomas Jahns <jahns@dkrz.de>
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>

/* read until end of line or end of file */
ssize_t
getline(char **linebuf, size_t *linebuf_size, FILE *fp)
{
  size_t len, old_len, buf_size;
  char *buf;
  ssize_t status = 0;

  if (feof(fp)) return -1;

  if (!linebuf_size || !linebuf)
    {
      errno = EINVAL;
      return -1;
    }
  if (!*linebuf)
    {
      buf_size = 80;
      buf = malloc(buf_size);
      if (!buf) return -1;
    }
  else
    {
      buf_size = *linebuf_size;
      buf = *linebuf;
    }
  /* initially buf_size is guaranteed < INT_MAX */
  if ((status = (fgets(buf, (int) buf_size, fp) != NULL)))
    {
      len = strlen(buf);
      if (buf[len - 1] != '\n')
        {
          do
            {
              size_t increment = 2 * buf_size - len <= INT_MAX ? 2 * buf_size - len : INT_MAX - (buf_size - len);
              void *temp = realloc(buf, buf_size + increment);
              if (!temp)
                {
                  status = -1;
                  break;
                }
              buf = temp;
              buf_size += increment;
              old_len = len;
            }
          while (fgets(buf + len, (int) (buf_size - len), fp) && (len += strlen(buf + len)) > old_len && buf[len - 1] != '\n');
        }
      status = len > SSIZE_MAX ? SSIZE_MAX : (ssize_t) len;
    }

  *linebuf_size = buf_size;
  *linebuf = buf;

  return status;
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
