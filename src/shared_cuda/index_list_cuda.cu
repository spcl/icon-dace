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

#include "index_list_cuda.h"

#ifndef __HIP__
// CUDA implementation

#include <cuda.h>
#include <cub/device/device_select.cuh>
#include <cub/iterator/counting_input_iterator.cuh>

#include <unordered_map>
#include <memory>

// Use stream-ordered allocs if we're capturing a graph
namespace {
    bool isStreamCapturing(gpuStream_t stream) {
        cudaStreamCaptureStatus captureStatus;
        cudaStreamIsCapturing(stream, &captureStatus);
        return captureStatus != cudaStreamCaptureStatusNone;
    }
}

class Storage {
public:
    virtual void  requestSize(size_t requestedSize) = 0;
    int* getNvalidPtr() {
        return reinterpret_cast<int*>(data);
    }
    char* getScratchPtr() {
        return data + alignment;
    }
    virtual ~Storage() = default;

protected:
    char* data = nullptr;
    static const int alignment = 512;
};
class AsyncStorage: public Storage {
public:
    AsyncStorage(gpuStream_t stream) :
      stream(stream) {    }

    void requestSize(size_t requestedSize) override final {
        if (data != nullptr) {
            cudaFreeAsync(data, stream);
        }
        cudaMallocAsync(&data, alignment+requestedSize, stream);
    }
    ~AsyncStorage() override {
        cudaFreeAsync(data, stream);
    }

private:
    gpuStream_t stream;
};

class SyncStorage : public Storage {
public:
    void requestSize(size_t requestedSize) override final {
        if (curSize < requestedSize+alignment) {
            cudaFree(data);
            cudaMalloc(&data, requestedSize+alignment);
            curSize = requestedSize+alignment;
        }
    }
    ~SyncStorage() override {
        cudaFree(data);
    }

private:
    size_t curSize = 0;
};

std::unordered_map<gpuStream_t, std::shared_ptr<SyncStorage>> syncStorageMap;


template<typename T>
struct ZeroCmp
{
    const T* conditions;
    const int startid;

    ZeroCmp(const int startid, const T* conditions) :
        startid(startid), conditions(conditions)
    { }

    __device__ __host__ __forceinline__
    bool operator() (const int &id)
    {
      return (conditions[ id - startid ] != 0);
    }
};

template <typename T>
static
void c_generate_index_list_gpu_generic_device(
            const T* dev_conditions,
            const int startid, const int endid,
            int* dev_indices,
            int* dev_nvalid, gpuStream_t stream)
{
    const int n = endid - startid + 1;

    // Argument is the offset of the first element
    cub::CountingInputIterator<int> iterator(startid);

    // Determine temporary device storage requirements
    size_t storageRequirement;
    cub::DeviceSelect::Flagged(nullptr, storageRequirement,
            iterator, dev_conditions, dev_indices,
            dev_nvalid, n, stream);

    // Allocate temporary storage
    // Use async storage in case we're capturing a graph
    // otherwise the sync storage per-stream
    std::shared_ptr<Storage> storage;
    if (isStreamCapturing(stream)) {
        storage = std::make_shared<AsyncStorage>(stream);
    } else {
        if (syncStorageMap.find(stream) == syncStorageMap.end()) {
            syncStorageMap[stream] = std::make_shared<SyncStorage>();
        }
        storage = syncStorageMap[stream];
    }

    storage->requestSize(storageRequirement);
    if (dev_nvalid == nullptr) {
        dev_nvalid = storage->getNvalidPtr();
    }

    ZeroCmp<T> select(startid, dev_conditions);
    cub::DeviceSelect::If(
            storage->getScratchPtr(), storageRequirement,
            iterator, dev_indices,
            dev_nvalid, n,
            select, stream);
}


template <typename T>
static
void c_generate_index_list_gpu_batched_generic(
            const int batch_size,
            const T* dev_conditions, const int cond_stride,
            const int startid, const int endid,
            int* dev_indices, const int idx_stride,
            int* dev_nvalid, gpuStream_t stream)
{
    for (int i = 0; i < batch_size; i++)
        c_generate_index_list_gpu_generic_device(
                dev_conditions + cond_stride*i,
                startid, endid,
                dev_indices + idx_stride*i,
                dev_nvalid + i, stream);
}

template <typename T>
static
void c_generate_index_list_gpu_generic(
            const T* dev_conditions,
            const int startid, const int endid,
            int* dev_indices, int* ptr_nvalid,
            bool copy_to_host, gpuStream_t stream)
{
    int* local_dev_nvalid = nullptr;

    c_generate_index_list_gpu_generic_device(
            dev_conditions, startid, endid, dev_indices,
            copy_to_host ? local_dev_nvalid : ptr_nvalid, stream);

    if (copy_to_host) {
        cudaMemcpyAsync(ptr_nvalid, local_dev_nvalid, sizeof(int), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
    }
}

///
/// Exposed functions
/// 
/// Non-batched first
/// 
void c_generate_index_list_gpu_single(
            const void* dev_conditions,
            const int startid, const int endid,
            int* dev_indices, int* nvalid,
            int data_size, bool copy_to_host,
            gpuStream_t stream)
{
    switch (data_size) {
        case 1: 
            c_generate_index_list_gpu_generic(
                static_cast<const char*>(dev_conditions),
                startid, endid, dev_indices, nvalid, copy_to_host, stream);
            break;
        case 4: 
            c_generate_index_list_gpu_generic(
                static_cast<const int*> (dev_conditions),
                startid, endid, dev_indices, nvalid, copy_to_host, stream);
            break;
    }
}

/// 
/// And now batched
/// 
void c_generate_index_list_gpu_batched(
        const int batch_size,
        const void* dev_conditions, const int cond_stride,
        const int startid, const int endid,
        int* dev_indices, const int idx_stride,
        int* dev_nvalid, int data_size,
        gpuStream_t stream)
{
    switch (data_size) {
        case 1: 
            c_generate_index_list_gpu_batched_generic(
                    batch_size,
                    static_cast<const char*>(dev_conditions),
                    cond_stride,
                    startid, endid,
                    dev_indices, idx_stride,
                    dev_nvalid, stream);
            break;

        case 4: 
            c_generate_index_list_gpu_batched_generic(
                    batch_size,
                    static_cast<const int*> (dev_conditions),
                    cond_stride,
                    startid, endid,
                    dev_indices, idx_stride,
                    dev_nvalid, stream);
            break;
    }
}

#else
// HIP implementation

#include <hip/hip_runtime.h>
#include <hipcub/device/device_select.hpp>
#include <hipcub/iterator/counting_input_iterator.hpp>

#include <unordered_map>
#include <memory>


class Storage {
public:
    virtual void  requestSize(size_t requestedSize) = 0;
    int* getNvalidPtr() {
        return reinterpret_cast<int*>(data);
    }
    char* getScratchPtr() {
        return data + alignment;
    }
    virtual ~Storage() = default;

protected:
    char* data = nullptr;
    static const int alignment = 512;
};
class SyncStorage : public Storage {
public:
    void requestSize(size_t requestedSize) override final {
        if (curSize < requestedSize+alignment) {
            hipFree(data);
            hipMalloc(&data, requestedSize+alignment);
            curSize = requestedSize+alignment;
        }
    }
    ~SyncStorage() override {
        hipFree(data);
    }

private:
    size_t curSize = 0;
};



SyncStorage storage;


template<typename T>
struct ZeroCmp
{
    const T* conditions;
    const int startid;

    ZeroCmp(const int startid, const T* conditions) :
        startid(startid), conditions(conditions)
    { }

    __device__ __host__ __forceinline__
    bool operator() (const int &id)
    {
      return (conditions[ id - startid ] != 0);
    }
};

template <typename T>
static
void c_generate_index_list_gpu_generic_device(
            const T* dev_conditions,
            const int startid, const int endid,
            int* dev_indices,
            int* dev_nvalid, gpuStream_t stream)
{
    const int n = endid - startid + 1;

    // Argument is the offset of the first element
    hipcub::CountingInputIterator<int> iterator(startid);

    // Determine temporary device storage requirements
    size_t storageRequirement;
    hipcub::DeviceSelect::Flagged(nullptr, storageRequirement,
            iterator, dev_conditions, dev_indices,
            dev_nvalid, n, 0);

    // Allocate temporary storage
    storage.requestSize(storageRequirement);
    if (dev_nvalid == nullptr) {
        dev_nvalid = storage.getNvalidPtr();
    }

    ZeroCmp<T> select(startid, dev_conditions);
    hipcub::DeviceSelect::If(
            storage.getScratchPtr(), storageRequirement,
            iterator, dev_indices,
            dev_nvalid, n,
            select, 0);
}


template <typename T>
static
void c_generate_index_list_gpu_batched_generic(
            const int batch_size,
            const T* dev_conditions, const int cond_stride,
            const int startid, const int endid,
            int* dev_indices, const int idx_stride,
            int* dev_nvalid, gpuStream_t stream)
{
    for (int i = 0; i < batch_size; i++)
        c_generate_index_list_gpu_generic_device(
                dev_conditions + cond_stride*i,
                startid, endid,
                dev_indices + idx_stride*i,
                dev_nvalid + i, 0);
}

template <typename T>
static
void c_generate_index_list_gpu_generic(
            const T* dev_conditions,
            const int startid, const int endid,
            int* dev_indices, int* ptr_nvalid,
            bool copy_to_host, gpuStream_t stream)
{
    int* local_dev_nvalid = nullptr;

    c_generate_index_list_gpu_generic_device(
            dev_conditions, startid, endid, dev_indices,
            copy_to_host ? local_dev_nvalid : ptr_nvalid, 0);

    if (copy_to_host) {
        hipMemcpyAsync(ptr_nvalid, local_dev_nvalid, sizeof(int), hipMemcpyDeviceToHost, 0);
        hipStreamSynchronize(0);
    }
}

///
/// Exposed functions
/// 
/// Non-batched first
/// 
void c_generate_index_list_gpu_single(
            const void* dev_conditions,
            const int startid, const int endid,
            int* dev_indices, int* nvalid,
            int data_size, bool copy_to_host,
            gpuStream_t stream)
{
    switch (data_size) {
        case 1: 
            c_generate_index_list_gpu_generic(
                static_cast<const char*>(dev_conditions),
                startid, endid, dev_indices, nvalid, copy_to_host, 0);
            break;
        case 4: 
            c_generate_index_list_gpu_generic(
                static_cast<const int*> (dev_conditions),
                startid, endid, dev_indices, nvalid, copy_to_host, 0);
            break;
    }
}

/// 
/// And now batched
/// 
void c_generate_index_list_gpu_batched(
        const int batch_size,
        const void* dev_conditions, const int cond_stride,
        const int startid, const int endid,
        int* dev_indices, const int idx_stride,
        int* dev_nvalid, int data_size,
        gpuStream_t stream)
{
    switch (data_size) {
        case 1: 
            c_generate_index_list_gpu_batched_generic(
                    batch_size,
                    static_cast<const char*>(dev_conditions),
                    cond_stride,
                    startid, endid,
                    dev_indices, idx_stride,
                    dev_nvalid, 0);
            break;

        case 4: 
            c_generate_index_list_gpu_batched_generic(
                    batch_size,
                    static_cast<const int*> (dev_conditions),
                    cond_stride,
                    startid, endid,
                    dev_indices, idx_stride,
                    dev_nvalid, 0);
            break;
    }
}

void initHIP(int deviceNum){
	hipSetDevice(deviceNum);
}
#endif 
