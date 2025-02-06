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

#ifndef UTIL_SYSINFO_H
#define UTIL_SYSINFO_H

#ifdef __cplusplus
extern "C" {
#endif

void util_user_name(char *name, int *actual_len);
void util_os_system(char *name, int *actual_len);
void util_node_name(char *name, int *actual_len);
void util_get_maxrss(int *maxrss);
void util_compiler_release(char *release_str, int *rstr_len);
void util_c_getpid(long int *pid);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_SYSINFO_H