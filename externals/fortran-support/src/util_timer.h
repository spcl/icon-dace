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

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#ifdef __cplusplus
extern "C" {
#endif

int util_cputime(double *user_time, double *system_time);
double util_walltime();
double util_gettimeofday();
void util_init_real_time();
void util_get_real_time_size(int *rt_size);
void util_read_real_time(void *it);
void util_diff_real_time(void *it1, void *it2, double *t);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_TIMER_H