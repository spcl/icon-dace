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

#pragma once

// Fortan interface to the following functions is
// implemented in ../src/shared/mo_index_list.f90

#ifdef __HIP__
#include <iostream>
#include <hip/hip_runtime.h>
using gpuStream_t = hipStream_t;
#else
// CUDA
using gpuStream_t = cudaStream_t;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

	void c_generate_index_list_gpu_single(
			const void* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int* nvalid, 
			int data_size, bool copy_to_host,
			gpuStream_t stream);

	void c_generate_index_list_gpu_batched(
			const int batch_size,
			const void* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, int data_size,
			gpuStream_t stream);

#ifdef __cplusplus
}
#endif
