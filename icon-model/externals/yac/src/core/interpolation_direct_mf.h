// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_DIRECT_MF_H
#define INTERPOLATION_DIRECT_MF_H

#include "interpolation_internal.h"

struct yac_interpolation_type * yac_interpolation_direct_mf_new(
  size_t collection_size, Xt_redist * redists, size_t num_src_fields);

#endif // INTERPOLATION_DIRECT_MF_H
