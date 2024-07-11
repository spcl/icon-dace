/**
 * @file xt_gpu.h
 * @brief routines for using GPU devices
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
 *
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

#ifndef XT_GPU_H
#define XT_GPU_H

#ifdef XT_GPU_INSTR_ENABLE

#define XT_GPU_INSTR_PUSH(arg) xt_gpu_instr_push(#arg)
#define XT_GPU_INSTR_POP xt_gpu_instr_pop()

#else

#define XT_GPU_INSTR_PUSH(arg)
#define XT_GPU_INSTR_POP

#endif

#include <stddef.h>

#include "core/ppm_visibility.h"

enum xt_memtype {
  XT_MEMTYPE_HOST,
  XT_MEMTYPE_DEVICE,
  XT_MEMTYPE_COUNT,
};

struct xt_gpu_vtable {
  void* (*Malloc)(size_t alloc_size, enum xt_memtype memtype);
  void (*Free)(void * ptr, enum xt_memtype memtype);
  void (*Memcpy)(
    void * dst, void const * src, size_t buffer_size,
    enum xt_memtype dst_memtype, enum xt_memtype src_memtype);
  enum xt_memtype (*Get_memtype)(const void *ptr);
  int (*Instr_push)(char const * name);
  int (*Instr_pop)(void);
};

/**
 * initialises xt_gpu
 * @remark if yaxt was compiled with CUDA support, xt_gpu will try to use it
 */
PPM_DSO_INTERNAL void
xt_gpu_init(void);

/**
 * allocates memory of the specified type
 * @param[in] alloc_size number of bytes to be allocated
 * @param[in] memtype    type of memory to be allocated
 * @return allocated memory of specified size and memory type\n
           NULL, if allocation failed
 */
PPM_DSO_INTERNAL void *
xt_gpu_malloc(size_t alloc_size, enum xt_memtype memtype);

/**
 * frees memory that was previously allocated using \ref xt_gpu_malloc
 * @param[in] ptr     pointer to be freed
 * @param[in] memtype type memory associated to ptr
 */
PPM_DSO_INTERNAL void
xt_gpu_free(void * ptr, enum xt_memtype memtype);

/**
 * copies memory from src to dst
 * @param[in] dst         pointer to destination memory
 * @param[in] src         pointer to source memory
 * @param[in] buffer_size number of bytes to be copied
 * @param[in] dst_memtype type of destination memory
 * @param[in] src_memtype type of source memory
 */
PPM_DSO_INTERNAL void xt_gpu_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_memtype dst_memtype, enum xt_memtype src_memtype);

/**
 * determines the type of memory associated with the provided pointer
 * @param[in] ptr pointer to be checked
 * @return type of memory associated to ptr
 */
enum xt_memtype xt_gpu_get_memtype(const void *ptr);

PPM_DSO_INTERNAL int xt_gpu_instr_push(char const * name);
PPM_DSO_INTERNAL int xt_gpu_instr_pop(void);

#endif // XT_GPU_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
