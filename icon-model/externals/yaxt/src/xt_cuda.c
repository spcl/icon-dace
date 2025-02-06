/**
 * @file xt_cuda.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "xt_cuda.h"

#include <cuda.h>
#include <dlfcn.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_gpu.h"

#ifdef XT_CUDA_NVTX_ENABLE
#include <nvToolsExt.h>
#endif

static const char filename[] = "xt_cuda.c";

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

#define DLSYM(name)                                         \
do {                                                        \
  void * ptr = dlsym(libcuda_handle, STRINGIFY(name));      \
  if (ptr != NULL) {                                        \
    *(void **)(&(func_ ## name)) = ptr;                     \
  } else {                                                  \
    load_successful = false;                                \
    fprintf(                                                \
      stderr, "ERROR:failed to load routine \"%s\" "        \
      "from %s\n", STRINGIFY(name), lib);                   \
  }                                                         \
} while(0)

#define CU_ERROR_CHECK(ret)                             \
do{                                                     \
    CUresult err = ret;                                 \
    if(err != CUDA_SUCCESS)                             \
    {                                                   \
      char const * err_string;                          \
      if (CUDA_SUCCESS !=                               \
            func_cuGetErrorString(err, &err_string))    \
        err_string = "undefined error";                 \
      fprintf(stderr, "Cuda driver error %d %s:: %s\n", \
              __LINE__, __func__, err_string);          \
      exit(EXIT_FAILURE);                               \
    }                                                   \
} while(0)

static CUresult (*func_cuGetErrorString)(CUresult error, const char** pStr);
static CUresult (*func_cuPointerGetAttribute)
  (void* data, CUpointer_attribute attribute, CUdeviceptr ptr);
static CUresult (*func_cuMemAlloc)(CUdeviceptr* dptr, size_t bytesize);
static CUresult (*func_cuMemFree)(CUdeviceptr dptr);
static CUresult (*func_cuMemcpyDtoD)(
  CUdeviceptr dstDevice, CUdeviceptr srcDevice, size_t ByteCount);
static CUresult (*func_cuMemcpyHtoD)(
  CUdeviceptr dstDevice, const void* srcHost, size_t ByteCount);
static CUresult (*func_cuMemcpyDtoH)(
  void* dstHost, CUdeviceptr srcDevice, size_t ByteCount);

static void * xt_cuda_malloc(
  size_t alloc_size, enum xt_memtype memtype) {

  switch (memtype) {
    default:
      Xt_abort(
        Xt_default_comm,
        "ERROR(xt_cuda_malloc): unsupported memory type",
        filename, __LINE__);
      return NULL;
    case (XT_MEMTYPE_HOST):
      return xmalloc(alloc_size);
    case (XT_MEMTYPE_DEVICE): {
      CUdeviceptr dptr;
      // cuMemAlloc fails for (alloc_size == 0)
      CU_ERROR_CHECK(func_cuMemAlloc(&dptr, MAX(alloc_size,1)));
      return (void*)dptr;
    }
  }
}

static void xt_cuda_free(void * ptr, enum xt_memtype memtype) {

  switch (memtype) {
    default:
      Xt_abort(
        Xt_default_comm,
        "ERROR(xt_cuda_free): unsupported memory type",
        filename, __LINE__);
      break;;
    case (XT_MEMTYPE_HOST):
      free(ptr);
      break;
    case (XT_MEMTYPE_DEVICE):
      CU_ERROR_CHECK(func_cuMemFree((CUdeviceptr)ptr));
      break;
  }
}

static void xt_cuda_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_memtype dst_memtype, enum xt_memtype src_memtype) {

  if (src_memtype == dst_memtype) {
    switch (src_memtype) {
      default:
        Xt_abort(
          Xt_default_comm,
          "ERROR(xt_cuda_memcpy): unsupported memory type",
          filename, __LINE__);
        break;
      case (XT_MEMTYPE_HOST):
        memcpy(dst, src, buffer_size);
        break;
      case (XT_MEMTYPE_DEVICE):
        CU_ERROR_CHECK(
          func_cuMemcpyDtoD(
            (CUdeviceptr)dst, (CUdeviceptr)src, buffer_size));
        break;
    }
  } else {
    switch (src_memtype) {
      default:
        Xt_abort(
          Xt_default_comm,
          "ERROR(xt_cuda_memcpy): unsupported source memory type",
          filename, __LINE__);
        break;
      case (XT_MEMTYPE_HOST):
        if (dst_memtype != XT_MEMTYPE_DEVICE)
          Xt_abort(
            Xt_default_comm,
            "ERROR(xt_cuda_memcpy): unsupported destination memory type",
            filename, __LINE__);
        CU_ERROR_CHECK(
          func_cuMemcpyHtoD((CUdeviceptr)dst, src, buffer_size));
        break;
      case (XT_MEMTYPE_DEVICE):
        if (dst_memtype != XT_MEMTYPE_HOST)
          Xt_abort(
            Xt_default_comm,
            "ERROR(xt_cuda_memcpy): unsupported destination memory type",
            filename, __LINE__);
        CU_ERROR_CHECK(
          func_cuMemcpyDtoH(dst, (CUdeviceptr)src, buffer_size));
        break;
    }
  }
}

static enum xt_memtype xt_cuda_get_memtype(const void *ptr) {

  XT_GPU_INSTR_PUSH(xt_cuda_get_memtype);

  CUmemorytype memorytype;
  CUresult ret =
    func_cuPointerGetAttribute(
      &memorytype, CU_POINTER_ATTRIBUTE_MEMORY_TYPE, (CUdeviceptr)ptr);

  XT_GPU_INSTR_POP;

  // in case the call to cuPointerGetAttribute failed, we just assume that
  // the pointer is a host pointer
  return
    ((ret == CUDA_SUCCESS) &&
     (memorytype == CU_MEMORYTYPE_DEVICE))?
       (XT_MEMTYPE_DEVICE):(XT_MEMTYPE_HOST);
}

#ifndef XT_CUDA_NVTX_ENABLE
static int dummy_instr_push(char const * XT_UNUSED(name)) {return 0;}
static int dummy_instr_pop() {return 0;}
#endif

static int load_cuda_library(void) {

  static int first_call = 1;
  static bool load_successful = false;

  if (!first_call) return load_successful;
  first_call = 0;

  static const char lib[] =
#ifdef __linux__
    "libcuda.so.1"
#elif defined __APPLE__ && defined __MACH__
    "libcuda.dylib"
#else
#warning "unsupported system, but trying libcuda.so.1"
    "libcuda.so.1"
#endif
    ;
  void *libcuda_handle = dlopen(lib, RTLD_NOW);

  load_successful = libcuda_handle != NULL;

  if (load_successful) {

    // load required CUDA library functions
    DLSYM(cuGetErrorString);
    DLSYM(cuPointerGetAttribute);
    DLSYM(cuMemAlloc);
    DLSYM(cuMemFree);
    DLSYM(cuMemcpyDtoD);
    DLSYM(cuMemcpyHtoD);
    DLSYM(cuMemcpyDtoH);

    dlclose(libcuda_handle);

  } else {

    int print_warning = 1;

    // check whether we are supposed to write a warning or not
    char const * cuda_warning_env =
      getenv("XT_CUDA_WARN_ON_MISSING_LIBCUDA");

    if (cuda_warning_env) {
      if (!strcmp(cuda_warning_env, "0")) print_warning = 0;
      else if (!strcmp(cuda_warning_env, "1")) print_warning = 1;
      else {
        Xt_abort(
          Xt_default_comm,
          "invalid value of XT_CUDA_WARN_ON_MISSING_LIBCUDA "
          "environment variable (has to be \"0\" or \"1\")",
          filename, __LINE__);
      }
    }

    if (print_warning) {
      // warn user about failed loading of CUDA library
      fputs(
        "-----------------------------------------------------------------------\n"
        "WARNING: yaxt was compiled with CUDA-support, but the library could not\n"
        "         be loaded. CUDA-support will be deactivated. Try setting\n"
        "         LD_LIBRARY_PATH to the location of libcuda.so.1 or set RPATH\n"
        "         accordingly.\n"
        "         To suppress this message set the\n"
        "         XT_CUDA_WARN_ON_MISSING_LIBCUDA environment variable to \"0\"\n"
        "-----------------------------------------------------------------------\n",
        stderr);
    }
  }

  return load_successful;
}

static struct xt_gpu_vtable const cuda_vtable = {
  .Malloc      = xt_cuda_malloc,
  .Free        = xt_cuda_free,
  .Memcpy      = xt_cuda_memcpy,
  .Get_memtype = xt_cuda_get_memtype,
#ifdef XT_CUDA_NVTX_ENABLE
  .Instr_push =  nvtxRangePush,
  .Instr_pop =   nvtxRangePop,
#else
  .Instr_push =  dummy_instr_push,
  .Instr_pop =   dummy_instr_pop,
#endif
};

const struct xt_gpu_vtable *xt_cuda_init(void) {

  return load_cuda_library() ? &cuda_vtable : (struct xt_gpu_vtable *)NULL;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
