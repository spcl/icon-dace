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
// implemented in ../src/shared/mo_util_system.f90

#include <stdlib.h>

#ifdef __XT3__
static int bsize = 65536;
static char *obuf, *ebuf;

void util_base_iobuf(void) {
    int stat1, stat2;

    obuf = malloc(bsize);
    ebuf = malloc(bsize);

    stat1 = setvbuf(stdout, obuf, _IOFBF, bsize);
    stat2 = setvbuf(stderr, ebuf, _IOFBF, bsize);
}
#endif

void util_exit(int exit_no) { exit(exit_no); }

void util_abort(void) { exit(1); }

int util_system(char *s) {
#ifdef __XT3__
    return (-1);
#else
    return (system(s));
#endif
}
