/**
 * @file xt_core.h
 *
 * \brief base definitions header file
 *
 * contains types used throughout the library and general initialization
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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

#ifndef XT_CORE_H
#define XT_CORE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>
#include <limits.h>

#include <mpi.h>

/** \example test_initialized_finalized.c
 */
/** \example test_initialized_finalized_f.f90
 */

/**
 * Distributed elements are enumerated with numbers of this signed
 * integral type.
 */
typedef XT_INT Xt_int;
#define Xt_int_dt (XT_INT_MPIDT)
typedef unsigned XT_INT Xt_uint;
#define Xt_uint_dt (XT_UINT_MPIDT)
typedef uint64_t Xt_uid;
#define XT_INT_MAX CONF_XT_INT_MAX
#define XT_INT_MIN CONF_XT_INT_MIN
/**
 * Use this in calls to printf, scanf etc. for data of type Xt_int
 */
#define XT_INT_FMT CONF_XT_INT_FMT

typedef struct Xt_idxlist_ *Xt_idxlist;
typedef struct Xt_xmap_ *Xt_xmap;
typedef struct Xt_redist_ *Xt_redist;

/**
 * represents \a size number of positions starting at \a start, where
 * the positions are assumed to decrease if \a size < 0 and increase otherwise.
 * Meaning the range of positions starting with
 * start up to start + size - 1 i.e. [start,start+size) if size is positive and
 * start down to start + size + 1 i.e. (start+size,start] if size is negative
 */
struct Xt_pos_ext
{
  int start, size;
};

/**
 * represents \a size number of positions beginning with \a start, where
 * the positions are (start + i * stride) for i in [0,size)
 */
struct Xt_offset_ext
{
  int start, size, stride;
};

/**
 * represents \a size number of offsets beginning with \a start, where
 * the offsets are \$start + i * stride for i in [0,size)\$
 * In contrast to Xt_offset_ext which is meant to hold the start as
 * number of elements, this struct holds start as MPI_Aint which is
 * meant to hold memory offsets, typically counted in bytes.
 */
struct Xt_aoffset_ext
{
  MPI_Aint start;
  int size;
  MPI_Aint stride;
};

/**
 * initialize library
 * @param[in] default_comm communicator to use for collective aborts
 */
void
xt_initialize(MPI_Comm default_comm);

/**
 * finalize library
 * @note this call only deallocates resources allocated via
 * xt_initialize and not objects explicitly created by xt_*_new calls.
 */
void
xt_finalize(void);

/**
 * @return 1 if xt_initialize has been called, zero otherwise
 */
int
xt_initialized(void);

/**
 * @return 1 if xt_finalized has been called, zero otherwise
 */
int
xt_finalized(void);


#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
