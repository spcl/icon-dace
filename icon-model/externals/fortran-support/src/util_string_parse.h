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

#ifndef UTIL_STRING_PARSE_H
#define UTIL_STRING_PARSE_H

#ifdef __cplusplus
extern "C" {
#endif

void do_parse_intlist(const char *in_parse_line, const int nvalues,
                      const int nlev_val, int *out_values, int *ierr);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_STRING_PARSE_H
