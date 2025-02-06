#ifndef GET_NUM_MISSVALS_H
#define GET_NUM_MISSVALS_H

#include <stddef.h>

size_t get_num_missvalsSP(size_t size, float *data, float missval);
size_t get_num_missvalsDP(size_t size, double *data, double missval);
size_t get_cplx_num_missvalsSP(size_t size, float *data, float missval);
size_t get_cplx_num_missvalsDP(size_t size, double *data, double missval);

#endif
