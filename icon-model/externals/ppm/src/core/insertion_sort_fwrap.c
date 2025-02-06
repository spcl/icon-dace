/**
 * @file insertion_sort_fwrap.c
 * @brief perform insertion sort on array
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
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
#  include "config.h"
#endif

#include <stdlib.h>
#include "core/insertion_sort.h"
#include "cfortran.h"

static void
PPM_insertion_sort_wrapper(void *a, int n, int es,
                           void *data, PPM_CompareWithData cmp)
{
  PPM_insertion_sort(a, (size_t)n, (size_t)es, data, cmp);
}


#undef ROUTINE_5
#define ROUTINE_5 (PPM_CompareWithData)(void (*)(void))
FCALLSCSUB5(PPM_insertion_sort_wrapper, PPM_INSERTION_SORT,
            ppm_insertion_sort,
            PVOID, INT, INT, PVOID, ROUTINE)

static void
PPM_sorted_insertion_wrapper(void *a, int n, int es,
                             void *data, PPM_CompareWithData cmp)
{
  PPM_sorted_insertion(a, (size_t)n, (size_t)es, data, cmp);
}


FCALLSCSUB5(PPM_sorted_insertion_wrapper, PPM_SORTED_INSERTION,
            ppm_sorted_insertion,
            PVOID, INT, INT, PVOID, ROUTINE)

static void
PPM_insertion_sort_once_wrapper(void *a, int n, int es, int inversion,
                                void *data, PPM_CompareWithData cmp)
{
  PPM_insertion_sort_once(a, (size_t)n, (size_t)es, (size_t)(inversion - 1),
                          data, cmp);
}

#undef ROUTINE_6
#define ROUTINE_6 (PPM_CompareWithData)(void (*)(void))
FCALLSCSUB6(PPM_insertion_sort_once_wrapper, PPM_INSERTION_SORT_ONCE,
            ppm_insertion_sort_once,
            PVOID, INT, INT, INT, PVOID, ROUTINE)

/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
