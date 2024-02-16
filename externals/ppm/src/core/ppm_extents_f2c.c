/**
 * @file ppm_extents_f2c.c
 * @brief Makes extents C API functions available to Fortran
 *
 * @copyright Copyright  (C)  2017  Thomas Jahns <jahns@dkrz.de>
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
#include "config.h"
#endif

#include "core/ppm_visibility.h"
#define FCALLSC_QUALIFIER PPM_DSO_API_EXPORT
#include "cfortran.h"
#include "core/ppm_extents.h"

FCALLSCSUB2(PPM_sprint_extent,
            PPM_SPRINT_EXTENT, ppm_sprint_extent,
            PSTRING, PVOID);

FCALLSCSUB2(PPM_sprint_extent64,
            PPM_SPRINT_EXTENT_I8, ppm_sprint_extent_i8,
            PSTRING, PVOID);

FCALLSCSUB2(PPM_sprint_iinterval,
            PPM_SPRINT_IINTERVAL, ppm_sprint_iinterval,
            PSTRING, PVOID);

FCALLSCSUB2(PPM_sprint_iinterval64,
            PPM_SPRINT_IINTERVAL_I8, ppm_sprint_iinterval_i8,
            PSTRING, PVOID);

FCALLSCSUB2(PPM_sprint_iinterval_sp,
            PPM_SPRINT_IINTERVAL_SP, ppm_sprint_iinterval_sp,
            PSTRING, PVOID);

FCALLSCSUB2(PPM_sprint_iinterval_dp,
            PPM_SPRINT_IINTERVAL_DP, ppm_sprint_iinterval_dp,
            PSTRING, PVOID);

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
