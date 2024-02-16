/**
 * @file interp_method_nnn.h
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

#ifndef INTERP_METHOD_NNN_H
#define INTERP_METHOD_NNN_H

#include "interp_method.h"

/** \example test_interp_method_nnn_parallel.c
 * A test for the parallel nearest neighbour interpolation method.
 */

/** \example test_interp_method_nnn_parallel2.c
 * A test for the parallel nearest neighbour interpolation method.
 */

/** \example test_interp_method_rbf_parallel.c
 * A test for the parallel radial basis function interpolation method.
 */


#define YAC_DEFAULT_RBF_SCALE (1.487973e+01)
#define YAC_DEFAULT_GAUSS_SCALE (0.1)

enum yac_interp_nnn_weight_type {
  NNN_AVG   = 0, //!< average of n source points
  NNN_DIST  = 1, //!< distance weighted average of n source points
  NNN_GAUSS = 2, //!< distance with Gauss weights of n source points
  NNN_RBF   = 3, //!< radial basis functions
};

struct yac_nnn_config {
  enum yac_interp_nnn_weight_type type;
  size_t n;
  union {
    double rbf_scale;
    double gauss_scale;
  } data;
};

struct interp_method * yac_interp_method_nnn_new(struct yac_nnn_config config);

#endif // INTERP_METHOD_NNN_H
