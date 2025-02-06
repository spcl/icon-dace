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

#ifndef UTIL_HASH_H
#define UTIL_HASH_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t util_hashword(const void *key, size_t length, uint32_t initval);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_HASH_H
