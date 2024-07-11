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

#ifndef UTIL_SYSTEM_H
#define UTIL_SYSTEM_H

#ifdef __cplusplus
extern "C" {
#endif

void util_exit(int exit_no);
void util_abort(void);
int util_system(char *s);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_SYSTEM_H