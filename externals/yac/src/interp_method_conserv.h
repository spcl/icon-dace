/**
 * @file interp_method_conserv.h
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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
 */

#ifndef INTERP_METHOD_CONSERV_H
#define INTERP_METHOD_CONSERV_H

#include "interp_method.h"

/** \example test_interp_method_conserv_parallel.c
 * A test for the parallel conservative‚ interpolation method.
 */

/**
 * normalisation options for conservative interpolation\n
 * (see SCRIP user manual for a more detailed description of the options)
 */
enum yac_interp_method_conserv_normalisation {
   CONSERV_DESTAREA = 0, //!< interpolation values are normalised by the area of
                         //!< the respective target cell
                         //!< (this is the default option)
   CONSERV_FRACAREA = 1, //!< interpolation values will be normalised by the area
                         //!< of the respective target cell that is covered by
                         //!< non-masked target cells
};

/**
 * constructor for a interpolation method of type interp_method_conserv
 * @param[in] order            1st or 2nd order remapping
 * @param[in] enforced_conserv enforce conservation by correcting truncation errors
 * @param[in] partial_coverage if == 0, target cells that are not completely
 *                             covered by non-masked source cells will not be
 *                             interpolated (this is the default option)\n
 *                             if != 0, for target cells that are only
 *                             partially covered by non-mask source cells the
 *                             interpolation value will be determined by the
 *                             contributions of the respective source cells
 * @param[in] normalisation    specifies how to do the normalisation
 */
struct interp_method * yac_interp_method_conserv_new(
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation);

#endif // INTERP_METHOD_CONSERV_H
