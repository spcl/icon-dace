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

// Fortan interface to the following functions is
// implemented in ../src/shared/mo_util_stride.f90

#include <stddef.h>

void util_stride_1d(int *out, int elemsize, void *p1, void *p2) {
    ptrdiff_t d = (unsigned char *) p2 - (unsigned char *) p1;
    d /= elemsize;
    *out = (int) d;
}

void util_stride_2d(int *out, int elemsize, const void *p1, const void *p2,
                    const void *p3) {
    ptrdiff_t d = (unsigned char *) p2 - (unsigned char *) p1;
    out[0] = (int) (d / elemsize);
    d = (unsigned char *) p3 - (unsigned char *) p1;
    out[1] = (int) (d / elemsize);
}

size_t util_get_ptrdiff(const void *a, const void *b) {
    return (size_t) ((unsigned char *) a > (unsigned char *) b
                         ? (unsigned char *) a - (unsigned char *) b
                         : (unsigned char *) b - (unsigned char *) a);
}
