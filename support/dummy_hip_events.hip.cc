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

#include <hip/hip_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

hipError_t hipEventRecord(hipEvent_t, hipStream_t) {
  return hipSuccess;
}

hipError_t hipStreamWaitEvent(hipStream_t, hipEvent_t, unsigned int) {
  return hipSuccess;
}

#ifdef __cplusplus
}
#endif
