// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

#ifndef UTIL_STRIDE_H
#define UTIL_STRIDE_H

#ifdef __cplusplus
extern "C" {
#endif

void util_stride_1d(int *out, int elemsize, void *p1, void *p2);
void util_stride_2d(int *out, int elemsize, const void *p1, const void *p2,
                    const void *p3);
size_t util_get_ptrdiff(const void *a, const void *b);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_STRIDE_H
