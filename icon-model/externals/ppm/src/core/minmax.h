/**
 * @file minmax.h
 * @brief maximum and minimum macros for ScalES-PPM, adapted from genometools
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef MINMAX_H
#define MINMAX_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

/*
  This file defines macros for maximum and minimum computation,
  if they are not already defined.
*/

/*
 * Don't nest these macros, random chaos will ensue!
 * You have been warned.
 */
#ifndef MAX
#define MAX(a,b)                  \
  ({ __typeof (a) _a = (a);       \
    __typeof (b) _b = (b);        \
    _a > _b ? _a : _b; })
#endif

#ifndef MIN
#define MIN(a,b)                \
  ({ __typeof (a) _a = (a);       \
    __typeof (b) _b = (b);        \
    _a <= _b ? _a : _b; })
#endif

#ifndef MIN3
#define MIN3(a, b, c)                           \
  ({ __typeof(a) _a = (a);                      \
     __typeof(b) _b = (b);                      \
     __typeof(c) _c = (c);                      \
     (((_a)<(_b))?((_a)<(_c)?(_a):(_c)):((_b)<(_c)?(_b):(_c)));})
#endif

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
