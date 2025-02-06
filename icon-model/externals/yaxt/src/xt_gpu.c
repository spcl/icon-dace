/**
 * @file xt_gpu.c
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

#include <string.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_gpu.h"
#ifdef HAVE_CUDA
#  include "xt_cuda.h"
#endif

static void * default_malloc(size_t alloc_size, enum xt_memtype memtype);
static void default_free(void * ptr, enum xt_memtype memtype);
static void default_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_memtype dst_memtype, enum xt_memtype src_memtype);
static enum xt_memtype default_get_memtype(const void *ptr);
static int default_instr_push(char const *);
static int default_instr_pop(void);

static struct xt_gpu_vtable vtable = {
  .Malloc      = default_malloc,
  .Free        = default_free,
  .Memcpy      = default_memcpy,
  .Get_memtype = default_get_memtype,
  .Instr_push  = default_instr_push,
  .Instr_pop   = default_instr_pop,
};

void xt_gpu_init() {

#ifdef HAVE_CUDA
  struct xt_gpu_vtable const * cuda_vtable = xt_cuda_init();
  if (cuda_vtable != NULL) vtable = *cuda_vtable;
#endif
}

#define STR(s) #s
#define CHECK_FOR_HOST_MEM(type)                             \
  do {                                                       \
    if (type != XT_MEMTYPE_HOST) {                           \
      Xt_abort(                                              \
        Xt_default_comm,                                     \
        "ERROR(" STR(__func__) "): unsupported memory type", \
        __FILE__, __LINE__);                                 \
    }                                                        \
  } while (0);

static void * default_malloc(size_t alloc_size, enum xt_memtype memtype) {

  CHECK_FOR_HOST_MEM(memtype);

  return xmalloc(alloc_size);
}

static void default_free(void * ptr, enum xt_memtype memtype) {

  CHECK_FOR_HOST_MEM(memtype);

  free(ptr);
}

static void default_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_memtype dst_memtype, enum xt_memtype src_memtype) {

  CHECK_FOR_HOST_MEM(src_memtype);
  CHECK_FOR_HOST_MEM(dst_memtype);

  memcpy(dst, src, buffer_size);
}

static enum xt_memtype default_get_memtype(const void *ptr)
{
  (void)ptr;
  return XT_MEMTYPE_HOST;
}

static int default_instr_push(char const * XT_UNUSED(name)) {return 0;}
static int default_instr_pop(void) {return 0;}

void * xt_gpu_malloc(size_t alloc_size, enum xt_memtype memtype) {
  return vtable.Malloc(alloc_size, memtype);
}

void xt_gpu_free(void * ptr, enum xt_memtype memtype) {
  vtable.Free(ptr, memtype);
}

void xt_gpu_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_memtype dst_memtype, enum xt_memtype src_memtype) {
  vtable.Memcpy(dst, src, buffer_size, dst_memtype, src_memtype);
}

enum xt_memtype xt_gpu_get_memtype(const void *ptr) {
  return vtable.Get_memtype(ptr);
}

int xt_gpu_instr_push(char const * name) {
  return vtable.Instr_push(name);
}

int xt_gpu_instr_pop() {
  return vtable.Instr_pop();
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
