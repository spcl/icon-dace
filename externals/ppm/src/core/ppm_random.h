/**
 * @file ppm_random.h
 * @brief PRNG C interface
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: pseudo random number generator interface
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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

#ifndef PPM_RANDOM_H
#define PPM_RANDOM_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <limits.h>
#include <inttypes.h>
#include "core/ppm_extents.h"
#include "core/core.h"

unsigned int
PPM_ya_rand_init(MPI_Comm comm, int);
void PPM_ya_rand_finish(void);



#define IRAND_MIN (INT_MIN + 1)
#define IRAND_MAX INT_MAX

/**
 * PRNG function for type int
 * @return random number in range [0,IRAND_MAX]
 **/
int
PPM_irandp(void);

/**
 * PRNG function for type int
 * @return random number in range [IRAND_MIN,IRAND_MAX]
 **/
int
PPM_irand(void);

/**
 * PRNG function for type int uint32_t
 * @return random number in range [0,2^32-1]
 **/
uint32_t
PPM_ya_random(void);

/**
 * PRNG function for type int uint64_t
 * @return random number in range [0,2^64-1]
 **/
uint64_t
PPM_ya_random64(void);

/**
 * PRNG function for type int int64_t
 * @return random number in range [0,2^63-1]
 **/
int64_t
PPM_irandp8(void);

/**
 * PRNG function for type int int64_t
 * @return random number in range [-2^63+1,2^63-1]
 **/
int64_t
PPM_irand8(void);


/**
 * PRNG function for type int
 * @param range range in which to generate random number
 * @return random number in range [range.first,range.last]
 **/
int
PPM_irandr(struct PPM_iinterval range);

/**
 * PRNG function for type int64_t
 * @param range range in which to generate random number
 * @return random number in range [range.first,range.last]
 **/
int64_t
PPM_irandr8(struct PPM_iinterval64 range);

/**
 * PRNG function for array of type int
 * @param a pointer to array to fill with
 * random numbers in range [IRAND_MIN,IRAND_MAX]
 * @param n number of elements in a to fill
 **/
void
PPM_irand_a(int *a, size_t n);

/**
 * PRNG function for array of type int
 * @param a pointer to array to fill with random numbers in range [0,IRAND_MAX]
 * @param n number of elements in a to fill
 **/
void
PPM_irandp_a(int *a, size_t n);

/**
 * PRNG function for array of type int
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_irandr_a(int *a, size_t n, struct PPM_iinterval range);

/**
 * PRNG function for array of type int, must be called from all threads
 * of an OpenMP thread team
 * @param a pointer to array to fill with
 * random numbers in range [IRAND_MIN,IRAND_MAX]
 * @param n number of elements in a to fill
 **/
void
PPM_irand_mt_a(int *a, size_t n);

/**
 * PRNG function for array of type int, must be called from all threads
 * of an OpenMP thread team
 * @param a pointer to array to fill with
 * random numbers in range [0,IRAND_MAX]
 * @param n number of elements in a to fill
 **/
void
PPM_irandp_mt_a(int *a, size_t n);

/**
 * PRNG function for array of type int, must be called from all threads
 * of an OpenMP thread team
 * @param a pointer to array to fill with
 * random numbers in range [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_irandr_mt_a(int *a, size_t n, struct PPM_iinterval range);

/**
 * PRNG function for uniformly distributed double precision
 * floating-point quantities.
 * This routine is synonymous with PPM_ya_fsgrandom.
 * @return random number in range (-1.0,1.0)
 **/
double
PPM_drand(void);

/**
 * PRNG function for uniformly distributed double precision
 * floating-point quantities
 * This routine is synonymous with PPM_drand.
 * @return random number in range (-1.0,1.0)
 **/
double
PPM_ya_fsgrandom(void);

/**
 * PRNG function for uniformly distributed double precision
 * floating-point quantities.
 * This routine is synonymous with PPM_ya_frandom.
 * @return random number in range [0.0,1.0)
 **/
double
PPM_drandp(void);

/**
 * PRNG function for uniformly distributed double precision
 * floating-point quantities
 * This routine is synonymous with PPM_drandp.
 * @return random number in range [0.0,1.0)
 **/
double
PPM_ya_frandom(void);

/**
 * PRNG function for uniformly distributed double precision
 * floating-point quantities
 * This routine is synonymous with PPM_drandp.
 * @param range range in which to generate random numbers
 * @return random number in range [range.first,range.last]
 **/
double
PPM_drandr(struct PPM_iinterval_dp range);

/**
 * PRNG function for array of type double
 *
 * @param a pointer to array to fill with random numbers in range
 * (-1.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_drand_a(double *a, size_t n);

/**
 * PRNG function for array of type double, must be called from all threads
 * of an OpenMP thread team.
 *
 * @param a pointer to array to fill with random numbers in range
 * (-1.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_drand_mt_a(double *a, size_t n);

/**
 * PRNG function for array of type double
 *
 * @param a pointer to array to fill with random numbers in range
 * [0.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_drandp_a(double *a, size_t n);

/**
 * PRNG function for array of type double, must be called from all threads
 * of an OpenMP thread team
 *
 * @param a pointer to array to fill with random numbers in range
 * [0.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_drandp_mt_a(double *a, size_t n);

/**
 * PRNG function for array of type double
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_drandr_a(double *a, size_t n, struct PPM_iinterval_dp range);

/**
 * PRNG function for array of type double, must be called from all threads
 * of an OpenMP thread team
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_drandr_mt_a(double *a, size_t n, struct PPM_iinterval_dp range);

/**
 * PRNG function for uniformly distributed single precision
 * floating-point quantities.
 * This routine is synonymous with PPM_ya_fsgrandom.
 * @return random number in range (-1.0,1.0)
 **/
float
PPM_frand(void);

/**
 * PRNG function for uniformly distributed single precision
 * floating-point quantities
 * This routine is synonymous with PPM_frand.
 * @return random number in range (-1.0,1.0)
 **/
float
PPM_ya_fsgrandomf(void);

/**
 * PRNG function for uniformly distributed single precision
 * floating-point quantities.
 * This routine is synonymous with PPM_ya_frandom.
 * @return random number in range [0.0,1.0)
 **/
float
PPM_frandp(void);

/**
 * PRNG function for uniformly distributed single precision
 * floating-point quantities
 * This routine is synonymous with PPM_frandp.
 * @return random number in range [0.0,1.0)
 **/
float
PPM_ya_frandomf(void);

/**
 * PRNG function for uniformly distributed single precision
 * floating-point quantities
 * This routine is synonymous with PPM_frandp.
 * @param range range in which to generate random numbers
 * @return random number in range [range.first,range.last]
 **/
float
PPM_frandr(struct PPM_iinterval_sp range);

/**
 * PRNG function for array of type float
 *
 * @param a pointer to array to fill with random numbers in range
 * (-1.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_frand_a(float *a, size_t n);

/**
 * PRNG function for array of type float, must be called from all threads
 * of an OpenMP thread team.
 *
 * @param a pointer to array to fill with random numbers in range
 * (-1.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_frand_mt_a(float *a, size_t n);

/**
 * PRNG function for array of type float
 *
 * @param a pointer to array to fill with random numbers in range
 * [0.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_frandp_a(float *a, size_t n);

/**
 * PRNG function for array of type float, must be called from all threads
 * of an OpenMP thread team
 *
 * @param a pointer to array to fill with random numbers in range
 * [0.0,1.0)
 * @param n number of elements in a to fill
 **/
void
PPM_frandp_mt_a(float *a, size_t n);

/**
 * PRNG function for array of type float
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_frandr_a(float *a, size_t n, struct PPM_iinterval_sp range);

/**
 * PRNG function for array of type float, must be called from all threads
 * of an OpenMP thread team
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_frandr_mt_a(float *a, size_t n, struct PPM_iinterval_sp range);


/**
 * PRNG function for array of type int64_t
 * @param a pointer to array to fill with
 * random numbers in range [-2^63+1,2^63-1]
 * @param n number of elements in a to fill
 **/
void
PPM_irand8_a(int64_t *a, size_t n);

/**
 * PRNG function for array of type int64_t, must be called from all threads
 * of an OpenMP thread team
 * @param a pointer to array to fill with
 * random numbers in range [-2^63+1,2^63-1]
 * @param n number of elements in a to fill
 **/
void
PPM_irand8_mt_a(int64_t *a, size_t n);

/**
 * PRNG function for array of type int64_t
 * @param a pointer to array to fill with
 * random numbers in range [0,2^63-1]
 * @param n number of elements in a to fill
 **/
void
PPM_irandp8_a(int64_t *a, size_t n);

/**
 * PRNG function for array of type int64_t, must be called from all threads
 * of an OpenMP thread team
 * @param a pointer to array to fill with
 * random numbers in range [0,2^63-1]
 * @param n number of elements in a to fill
 **/
void
PPM_irandp8_mt_a(int64_t *a, size_t n);

/**
 * PRNG function for array of type int64_t
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_irandr8_a(int64_t *a, size_t n, struct PPM_iinterval64 range);

/**
 * PRNG function for array of type int64_t, must be called from all threads
 * of an OpenMP thread team
 *
 * @param a pointer to array to fill with random numbers in range
 * [range.first,range.last]
 * @param n number of elements in a to fill
 * @param range range in which to generate random numbers
 **/
void
PPM_irandr8_mt_a(int64_t *a, size_t n, struct PPM_iinterval64 range);


#endif

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
