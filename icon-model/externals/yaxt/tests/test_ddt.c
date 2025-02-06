/**
 * @file test_ddt.c
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
#include <config.h>
#endif

#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <complex.h>
#include <float.h>

#ifdef _OPENACC
#include <openacc.h>
#define STR(s) #s
#define my_Pragma(args) _Pragma(args)
#define PragmaACC(args) my_Pragma(STR(acc args))
#else
#define PragmaACC(args)
#endif

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "../src/xt_ddt.h"
#include "../src/xt_gpu.h"

#include "tests.h"
#include "ctest_common.h"

enum { padding_byte_value = -1 };

enum { ldbl_size = LDBL_MANT_DIG == 64 ? 10 : sizeof(long double) };

static int check_xt_ddt(
  MPI_Datatype mpi_ddt, void const * in_data, size_t in_data_size,
  size_t ref_pack_size);
static int check_xt_ddt_off(
  MPI_Datatype mpi_ddt, void const * in_data, int in_data_offset,
  size_t in_data_size, size_t ref_pack_size);

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);
  xt_gpu_init();

  { // elemental integer datatype
    MPI_Datatype mpi_ddt = MPI_INT;
    int const in_data[] = {1337};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR MPI_INT");
  }

  { // elemental double datatype
    MPI_Datatype mpi_ddt = MPI_DOUBLE;
    double const in_data[] = {1337.0};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR MPI_DOUBLE");
  }

  { // elemental float datatype
    MPI_Datatype mpi_ddt = MPI_FLOAT;
    float const in_data[] = {1337.0};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR MPI_FLOAT");
  }

#if (MPI_VERSION > 2 || ( MPI_VERSION == 2 && MPI_SUBVERSION >= 2)) && defined HAVE_LONG_DOUBLE__COMPLEX
#ifdef HAVE_DECL___BUILTIN_COMPLEX
#define Xt_complex_value(real,imag) __builtin_complex(real, imag)
#elif defined __clang__
#define Xt_complex_value(real,imag) (__extension__ (long double _Complex){(real), (imag)})
#else
#define Xt_complex_value(real,imag) (real + imag*I)
#endif
  if (MPI_C_LONG_DOUBLE_COMPLEX != MPI_DATATYPE_NULL)
  { // elemental long double complex datatype
    MPI_Datatype mpi_ddt = MPI_C_LONG_DOUBLE_COMPLEX;
    static const long double _Complex in_data[] = { Xt_complex_value(1337.0L, -1337.0L) };

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR MPI_C_LONG_DOUBLE_COMPLEX");
  }
#endif

  { // elemental float + int datatype
    MPI_Datatype mpi_ddt = MPI_FLOAT_INT;
    static const struct {
      float f;
      int i;
    } in_data[] = {{.f = 1337.0f, .i = 1337}};

    if (check_xt_ddt(
          mpi_ddt, in_data,
          sizeof(in_data), sizeof(in_data[0].f) + sizeof(in_data[0].i)))
      PUT_ERR("ERROR MPI_FLOAT_INT");
  }

  { // elemental double + int datatype
    MPI_Datatype mpi_ddt = MPI_DOUBLE_INT;
    static const struct {
      double d;
      int i;
    } in_data[] = {{.d = 1337.0, .i = 1337}};

    if (check_xt_ddt(
          mpi_ddt, in_data,
          sizeof(in_data), sizeof(in_data[0].d) + sizeof(in_data[0].i)))
      PUT_ERR("ERROR MPI_DOUBLE_INT");
  }

  if (MPI_LONG_DOUBLE_INT != MPI_DATATYPE_NULL)
  { // elemental long double + int datatype
    MPI_Datatype mpi_ddt = MPI_LONG_DOUBLE_INT;
    static const struct {
      long double ld;
      int i;
    } in_data[] = {{ .ld = 1337.0L, .i = 1337 }};
    if (check_xt_ddt(
          mpi_ddt, in_data,
          sizeof(in_data), sizeof(in_data[0].ld) + sizeof(in_data[0].i)))
      PUT_ERR("ERROR MPI_LONG_DOUBLE_INT");
  }

  { // MPI_Type_dup(MPI_INT)
    MPI_Datatype mpi_ddt;
    MPI_Type_dup(MPI_INT, &mpi_ddt);
    int const in_data[] = {1337};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR DUP(INT)");
  }

  { // MPI_Type_contiguous(7, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 7};
    MPI_Type_contiguous(COUNT, MPI_DOUBLE, &mpi_ddt);
    double in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), COUNT * sizeof(in_data[0])))
      PUT_ERR("ERROR CONT(7, DOUBLE)");
  }

  { // MPI_Dup(MPI_Type_contiguous(7, MPI_DOUBLE))
    MPI_Datatype mpi_ddt, cont_mpi_ddt;
    enum {COUNT = 7};
    MPI_Type_contiguous(COUNT, MPI_DOUBLE, &cont_mpi_ddt);
    MPI_Type_dup(cont_mpi_ddt, &mpi_ddt);
    MPI_Type_free(&cont_mpi_ddt);
    double in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), COUNT * sizeof(in_data[0])))
      PUT_ERR("ERROR DUP(CONT(7, DOUBLE))");
  }

#ifndef XT_CANNOT_SUPPORT_ZERO_SIZE_DT
  { // MPI_Type_contiguous(0, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 0};
    MPI_Type_contiguous(COUNT, MPI_DOUBLE, &mpi_ddt);
    double in_data[COUNT+1];
    for (size_t i = 0; i < COUNT+1; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), 0))
      PUT_ERR("ERROR CONT(0, DOUBLE)");
  }
#else
#warning "Skipping tests for defective MPI! Zero-size datatypes unsupported."
#endif                                          \

  { // MPI_Type_contiguous(1, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    MPI_Type_contiguous(COUNT, MPI_DOUBLE, &mpi_ddt);
    double in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), COUNT * sizeof(in_data[0])))
      PUT_ERR("ERROR CONT(1, DOUBLE)");
  }

  { // MPI_Type_contiguous(4, MPI_FLOAT_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    MPI_Type_contiguous(COUNT, MPI_FLOAT_INT, &mpi_ddt);
    struct {
      float f;
      int i;
    } in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) {
      in_data[i].f = (float)i;
      in_data[i].i = (int)i;
    }

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * (sizeof(in_data[0].f) + sizeof(in_data[0].i))))
      PUT_ERR("ERROR CONT(4, MPI_FLOAT_INT)");
  }

  { // MPI_Type_contiguous(4, MPI_DOUBLE_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    MPI_Type_contiguous(COUNT, MPI_DOUBLE_INT, &mpi_ddt);
    struct {
      double d;
      int i;
    } in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) {
      in_data[i].d = (double)i;
      in_data[i].i = (int)i;
    }

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * (sizeof(in_data[0].d) + sizeof(in_data[0].i))))
      PUT_ERR("ERROR CONT(4, MPI_DOUBLE_INT)");
  }

  if (MPI_LONG_DOUBLE_INT != MPI_DATATYPE_NULL)
  { // MPI_Type_contiguous(4, MPI_LONG_DOUBLE_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    MPI_Type_contiguous(COUNT, MPI_LONG_DOUBLE_INT, &mpi_ddt);
    struct {
      long double ld;
      int i;
    } in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) {
      in_data[i].ld = (long double)i;
      in_data[i].i = (int)i;
      // on x86 and x86_64 we need to fix the missing initialization
      // of the bytes following the 80 used bits of a long double
      memset((unsigned char *)(in_data+i)+ldbl_size, 0,
             sizeof (long double) - ldbl_size);
    }

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * (sizeof(in_data[0].ld) + sizeof(in_data[0].i))))
      PUT_ERR("ERROR CONT(4, MPI_LONG_DOUBLE_INT)");
  }

  { // MPI_Type_contiguous(3, MPI_Type_contiguous(3,MPI_FLOAT))
    enum {COUNT = 3};
    MPI_Datatype cont_mpi_ddt;
    MPI_Type_contiguous(COUNT, MPI_FLOAT, &cont_mpi_ddt);
    MPI_Datatype mpi_ddt;
    MPI_Type_contiguous(COUNT, cont_mpi_ddt, &mpi_ddt);
    MPI_Type_free(&cont_mpi_ddt);
    float in_data[COUNT][COUNT];
    for (size_t i = 0, k = 0; i < COUNT; ++i)
      for (size_t j = 0; j < COUNT; ++j, ++k)
        in_data[i][j] = (float)k;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * COUNT * sizeof(in_data[0][0])))
      PUT_ERR("ERROR CONT(3, CONT(3, DOUBLE))");
  }

  { // MPI_Type_vector(3, 5, 16, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 3, BLOCKLENGTH = 5, STRIDE = 16};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(3, 5, 16, FLOAT)");
  }

  { // MPI_Type_vector(1, 1, 16, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1, STRIDE = 16};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(1, 1, 16, FLOAT)");
  }

  if (MPI_LONG_DOUBLE != MPI_DATATYPE_NULL)
  { // MPI_Type_vector(3, 1, 16, MPI_LONG_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 3, BLOCKLENGTH = 1, STRIDE = 16};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_LONG_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    long double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) {
      in_data[i] = (long double)i;
      // On x86 and x86_64 we need to fix the missing initialization
      // of the bytes following the 80 bits actually used in a long
      // double. For alignment reasons, the sizeof (long double) is 12
      // or 16 on these platforms, even though only 10 are actually used.
      //
      // Also this needs to use the same byte value as used later to
      // initialize the buffer in the copy test because some MPI
      // implementations copy only the 10 bytes that matter and others
      // the whole 12 or 16, like xt_ddt does.
      memset((unsigned char *)(in_data+i)+ldbl_size, padding_byte_value,
             sizeof (long double) - ldbl_size);
    }

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(3, 1, 16, LONG_DOUBLE)");
  }

  { // MPI_Type_vector(1, 1, 1, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1, STRIDE = 1};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(1, 1, 1, FLOAT)");
  }

#ifndef XT_CANNOT_SUPPORT_ZERO_SIZE_DT
  { // MPI_Type_vector(0, 1, 1, MPI_CHAR)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 0, BLOCKLENGTH = 1, STRIDE = 1};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_CHAR, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    char in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (char)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH* sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(0, 1, 1, CHAR)");
  }
#endif

  { // MPI_Type_vector(3, 3, -3, MPI_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 3, BLOCKLENGTH = 3, STRIDE = -3};
    MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_INT, &mpi_ddt);
    enum {IN_DATA_COUNT = 9};
    int in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;

    if (check_xt_ddt_off(
          mpi_ddt, in_data, 6 * sizeof(int), sizeof(in_data),
          COUNT * BLOCKLENGTH* sizeof(in_data[0])))
      PUT_ERR("ERROR VECTOR(3, 3, -3, INT)");
  }

  { // MPI_Type_create_hvector(3, 5, 16*sizeof(int), MPI_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 3, BLOCKLENGTH = 5, STRIDE = 16*sizeof(int)};
    MPI_Type_create_hvector(COUNT, BLOCKLENGTH, STRIDE, MPI_INT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    int in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HVECTOR(3, 5, 16*sizeof(int), INT)");
  }

  { // MPI_Type_create_hvector(1, 1, 16*sizeof(int), MPI_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1, STRIDE = 16*sizeof(int)};
    MPI_Type_create_hvector(COUNT, BLOCKLENGTH, STRIDE, MPI_INT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    int in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HVECTOR(1, 1, 16*sizeof(int), INT)");
  }

  { // MPI_Type_create_hvector(1, 1, sizeof(int), MPI_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1, STRIDE = sizeof(int)};
    MPI_Type_create_hvector(COUNT, BLOCKLENGTH, STRIDE, MPI_INT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * STRIDE};
    int in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HVECTOR(1, 1, sizeof(int), INT)");
  }

  { // MPI_Type_indexed(4, {8,4,2,1}, {0,8,16,24}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    static const int array_of_blocklengths[COUNT] = {8,4,2,1};
    static const int array_of_displacements[COUNT] = {0,8,16,24};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 8};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 15 * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED(4, {8,4,2,1}, {0,8,16,24}, DOUBLE)");
  }

  { // MPI_Type_indexed(1, {1}, {2}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {1};
    static const int array_of_displacements[COUNT] = {2};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 3};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED(1, {1}, {2}, DOUBLE)");
  }

  { // MPI_Type_indexed(1, {1}, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {1};
    static const int array_of_displacements[COUNT] = {0};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 2};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED(1, {1}, {0}, DOUBLE)");
  }

#ifndef XT_CANNOT_SUPPORT_ZERO_SIZE_DT
  { // MPI_Type_indexed(0, {1}, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 0};
    static const int array_of_blocklengths[1] = {-1};
    static const int array_of_displacements[1] = {-1};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), 0))
      PUT_ERR("ERROR INDEXED(0, {1}, {0}, DOUBLE)");
  }

  { // MPI_Type_indexed(1, {0}, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {0};
    static const int array_of_displacements[COUNT] = {0};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), 0))
      PUT_ERR("ERROR INDEXED(1, {0}, {0}, DOUBLE)");
  }
#endif

  { // MPI_Type_create_hindexed(
    //   4, {8,4,2,1}, sizeof(float)*{0,8,16,24}, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    static const int array_of_blocklengths[COUNT] = {8,4,2,1};
    static const MPI_Aint array_of_displacements[COUNT] =
      {0*sizeof(float),8*sizeof(float),16*sizeof(float),24*sizeof(float)};
    MPI_Type_create_hindexed(
      COUNT, (int *)array_of_blocklengths, (MPI_Aint *)array_of_displacements,
      MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 8};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 15 * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED(4, {8,4,2,1}, sizeof(float)*{0,8,16,24}, FLOAT)");
  }

  { // MPI_Type_create_hindexed(1, {1}, sizeof(float)*{2}, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {1};
    static const MPI_Aint array_of_displacements[COUNT] = {2*sizeof(float)};
    MPI_Type_create_hindexed(
      COUNT, (int *)array_of_blocklengths, (MPI_Aint *)array_of_displacements,
      MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 4};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED(1, {1}, sizeof(float)*{2}, FLOAT)");
  }

  { // MPI_Type_create_hindexed(1, {1}, {0}, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {1};
    static const MPI_Aint array_of_displacements[COUNT] = {0};
    MPI_Type_create_hindexed(
      COUNT, (int *)array_of_blocklengths, (MPI_Aint *)array_of_displacements,
      MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = COUNT * 2};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED(1, {1}, {0}, FLOAT)");
  }

  { // MPI_Type_create_indexed_block(4, 3, {0,4,12,8}, MPI_INT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4, BLOCKLENGTH = 3};
    static const int array_of_displacements[COUNT] = {0,4,12,8};
    MPI_Type_create_indexed_block(
      COUNT, BLOCKLENGTH, (int *)array_of_displacements, MPI_INT, &mpi_ddt);
    enum {IN_DATA_COUNT = 16};
    int in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED_BLOCK(4, 3, {0,4,12,8}, INT)");
  }

  { // MPI_Type_create_indexed_block(1, 1, {2}, MPI_CHAR)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1};
    static const int array_of_displacements[COUNT] = {2};
    MPI_Type_create_indexed_block(
      COUNT, BLOCKLENGTH, (int *)array_of_displacements, MPI_CHAR, &mpi_ddt);
    enum {IN_DATA_COUNT = 3};
    char in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (char)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED_BLOCK(1, 1, {2}, CHAR)");
  }

  { // MPI_Type_create_indexed_block(1, 1, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1};
    static const int array_of_displacements[COUNT] = {0};
    MPI_Type_create_indexed_block(
      COUNT, BLOCKLENGTH, (int *)array_of_displacements, MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 2};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED_BLOCK(1, 1, {0}, DOUBLE)");
  }

#ifndef XT_CANNOT_SUPPORT_ZERO_SIZE_DT
  { // MPI_Type_create_indexed_block(0, 1, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 0, BLOCKLENGTH = 1};
    static const int array_of_displacements[1] = {0};
    MPI_Type_create_indexed_block(
      COUNT, BLOCKLENGTH, (int *)array_of_displacements, MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED_BLOCK(0, 1, {0}, DOUBLE)");
  }
#endif

  { // MPI_Type_create_indexed_block(1, 0, {0}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 0};
    static const int array_of_displacements[COUNT] = {0};
    MPI_Type_create_indexed_block(
      COUNT, BLOCKLENGTH, (int *)array_of_displacements, MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR INDEXED_BLOCK(1, 0, {0}, DOUBLE)");
  }

#if MPI_VERSION >= 3
  { // MPI_Type_create_hindexed_block(
    //   4, 3, sizeof(double)*{0,4,12,8}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 4, BLOCKLENGTH = 3};
    static const MPI_Aint array_of_displacements[COUNT] =
      {sizeof(double)*0,sizeof(double)*4,sizeof(double)*12,sizeof(double)*8};
    MPI_Type_create_hindexed_block(
      COUNT, BLOCKLENGTH, (MPI_Aint*)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 16};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED_BLOCK(4, 3, sizeof(double)*{0,4,12,8}, DOUBLE)");
  }

  { // MPI_Type_create_hindexed_block(
    //   1, 1, sizeof(double)*{2}, MPI_DOUBLE)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1};
    static const MPI_Aint array_of_displacements[COUNT] =
      {sizeof(double)*2};
    MPI_Type_create_hindexed_block(
      COUNT, BLOCKLENGTH, (MPI_Aint *)array_of_displacements,
      MPI_DOUBLE, &mpi_ddt);
    enum {IN_DATA_COUNT = 3};
    double in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED_BLOCK(1, 1, sizeof(double)*{2}, DOUBLE)");
  }

  { // MPI_Type_create_hindexed_block(
    //   1, 1, {0}, MPI_FLOAT)
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1, BLOCKLENGTH = 1};
    MPI_Aint const array_of_displacements[COUNT] = {0};
    MPI_Type_create_hindexed_block(
      COUNT, BLOCKLENGTH, array_of_displacements, MPI_FLOAT, &mpi_ddt);
    enum {IN_DATA_COUNT = 1};
    float in_data[IN_DATA_COUNT];
    for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          COUNT * BLOCKLENGTH * sizeof(in_data[0])))
      PUT_ERR("ERROR HINDEXED_BLOCK(1, 1, {0}, FLOAT)");
  }
#endif

  { // MPI_Type_create_struct(...)
    static const struct {
      char dummy1;
      int i[3];
      char dummy2[3];
      float f;
      char dummy3[5];
      double d[2];
    } in_data[] =
    {{.dummy1 = 'a',
      .i = {1, 2, 3},
      .dummy2 = {'b', 'c'},
      .f = 5.0,
      .dummy3 = {'d', 'e', 'f'},
      .d = {6.0, 7.0}}};
    MPI_Datatype mpi_ddt;
    enum {COUNT = 3};
    static const int array_of_blocklengths[COUNT] = {3, 1, 2};
    MPI_Aint base_address, i_address, f_address, d_address;
    MPI_Get_address((void *)in_data, &base_address);
    MPI_Get_address((void *)in_data[0].i, &i_address);
    MPI_Get_address((void *)&in_data[0].f, &f_address);
    MPI_Get_address((void *)in_data[0].d, &d_address);
    MPI_Aint array_of_displacements[COUNT] =
      {(MPI_Aint)(i_address - base_address),
       (MPI_Aint)(f_address - base_address),
       (MPI_Aint)(d_address - base_address)};
    static const MPI_Datatype array_of_types[COUNT] =
      {MPI_INT, MPI_FLOAT, MPI_DOUBLE};
    MPI_Type_create_struct(
      COUNT, (int *)array_of_blocklengths, (MPI_Aint *)array_of_displacements,
      (MPI_Datatype *)array_of_types, &mpi_ddt);

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          (3 * sizeof(int) + 1 * sizeof(float) + 2 * sizeof(double))))
      PUT_ERR("ERROR STRUCT(...)");
  }

  { // MPI_Type_create_struct(1, {1}, {0}, {MPI_INT})
    MPI_Datatype mpi_ddt;
    enum {COUNT = 1};
    static const int array_of_blocklengths[COUNT] = {1};
    static const MPI_Aint array_of_displacements[COUNT] = {0};
    static const MPI_Datatype array_of_types[COUNT] = {MPI_INT};
    int in_data[COUNT];
    for (size_t i = 0; i < COUNT; ++i) in_data[i] = (int)i;
    MPI_Type_create_struct(
      COUNT, (int *)array_of_blocklengths, (MPI_Aint *)array_of_displacements,
      (MPI_Datatype *)array_of_types, &mpi_ddt);

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR("ERROR STRUCT(1, {1}, {0}, {INT})");
  }

  { // MPI_Type_create_subarray(
    //   4, {8,4,3,7}, {2,3,1,4}, {3,1,0,2}, order, MPI_DOUBLE)
    int order_types[] = {MPI_ORDER_FORTRAN, MPI_ORDER_C};
    for (size_t order_idx = 0;
         order_idx < sizeof(order_types) / sizeof(order_types[0]);
         ++order_idx) {

      MPI_Datatype mpi_ddt;
      enum {NDIM = 4, DIM0 = 8, DIM1 = 4, DIM2 = 4, DIM3 = 7};
      static const int array_of_sizes[NDIM] = {DIM0, DIM1, DIM2, DIM3};
      static const int array_of_subsizes[NDIM] = {2, 3, 1, 4};
      static const int array_of_starts[NDIM] = {3,1,0,2};
      enum {IN_DATA_COUNT = DIM0 * DIM1 * DIM2 * DIM3};
      double in_data[IN_DATA_COUNT];
      for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (double)i;
      MPI_Type_create_subarray(
        NDIM, (int *)array_of_sizes, (int *)array_of_subsizes,
        (int *)array_of_starts, order_types[order_idx], MPI_DOUBLE, &mpi_ddt);

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data),
          2 * 3 * 1 * 4 * sizeof(in_data[0])))
      PUT_ERR(
        "ERROR SUBARRAY(4, {8,4,3,7}, {2,3,1,4}, {3,1,0,2}, order, DOUBLE)");
    }
  }

  { // MPI_Type_create_subarray(1, {8}, {4}, {2}, order, MPI_FLOAT)
    int order_types[] = {MPI_ORDER_FORTRAN, MPI_ORDER_C};
    for (size_t order_idx = 0;
         order_idx < sizeof(order_types) / sizeof(order_types[0]);
         ++order_idx) {

      MPI_Datatype mpi_ddt;
      enum {NDIM = 1, DIM0 = 8};
      static const int array_of_sizes[NDIM] = {DIM0};
      static const int array_of_subsizes[NDIM] = {4};
      static const int array_of_starts[NDIM] = {2};
      enum {IN_DATA_COUNT = DIM0};
      float in_data[IN_DATA_COUNT];
      for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (float)i;
      MPI_Type_create_subarray(
        NDIM, (int *)array_of_sizes, (int *)array_of_subsizes,
        (int *)array_of_starts, order_types[order_idx], MPI_FLOAT, &mpi_ddt);

      if (check_xt_ddt(
            mpi_ddt, in_data, sizeof(in_data), 4 * sizeof(in_data[0])))
        PUT_ERR(
          "ERROR SUBARRAY(1, {8}, {4}, {2}, order, FLOAT)");
    }
  }

  { // MPI_Type_create_subarray(1, {1}, {1}, {0}, order, MPI_INT)
    int order_types[] = {MPI_ORDER_FORTRAN, MPI_ORDER_C};
    for (size_t order_idx = 0;
         order_idx < sizeof(order_types) / sizeof(order_types[0]);
         ++order_idx) {

      MPI_Datatype mpi_ddt;
      enum {NDIM = 1, DIM0 = 1};
      static const int array_of_sizes[NDIM] = {DIM0};
      static const int array_of_subsizes[NDIM] = {1};
      static const int array_of_starts[NDIM] = {0};
      enum {IN_DATA_COUNT = DIM0};
      int in_data[IN_DATA_COUNT];
      for (size_t i = 0; i < IN_DATA_COUNT; ++i) in_data[i] = (int)i;
      MPI_Type_create_subarray(
        NDIM, (int *)array_of_sizes, (int *)array_of_subsizes,
        (int *)array_of_starts, order_types[order_idx], MPI_FLOAT, &mpi_ddt);

    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), 1 * sizeof(in_data[0])))
      PUT_ERR(
        "ERROR SUBARRAY(1, {1}, {1}, {0}, order, INT)");
    }
  }

  { // MPI_Type_create_resized(MPI_FLOAT, -3, 16)
    MPI_Datatype mpi_ddt;
    enum {LB = -3, EXTENT = 16};
    MPI_Type_create_resized(MPI_FLOAT, LB, EXTENT, &mpi_ddt);
    float in_data[] = {1337.0};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR RESIZED(MPI_FLOAT, -3, 16)");
  }

  { // MPI_Type_create_resized(MPI_INT, -4, 0)
    MPI_Datatype mpi_ddt;
    enum {LB = -4, EXTENT = 0};
    MPI_Type_create_resized(MPI_INT, LB, EXTENT, &mpi_ddt);
    int in_data[] = {1337};

    if (check_xt_ddt(mpi_ddt, in_data, sizeof(in_data), sizeof(in_data[0])))
      PUT_ERR("ERROR RESIZED(MPI_INT, -4, 0)");
  }

  { // INDEXED(
    //   4, {1,2,3,4}, {3,7,9,12},
    //   STRUCT(
    //     3,{1,1,1,1},{...},
    //     {CONT(0,DUP(MPI_INT)),
    //      CONT(3,MPI_INT16_T),
    //      VECTOR(3,3,9,MPI_DOUBLE),
    //      CONT(0,MPI_CHAR)}))
    // Unfortunately: on buggy OpenMPI versions, contiguous needs to
    //   have at least count 1, this seems to be fixed since OpenMPI
    //   4.0.2 and apparently wasn't present in very old releases.
#if defined OMPI_MAJOR_VERSION \
  && ( OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 8               \
       || OMPI_MAJOR_VERSION == 2                                       \
       || OMPI_MAJOR_VERSION == 3                                       \
       || OMPI_MAJOR_VERSION == 4 && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION == 1 )
    enum { minimal_cont_count = 1, };
#else
    enum { minimal_cont_count = 0, };
#endif
    struct {
      char c1[3];
      int i[2];
      int16_t i16[3];
      double d[9][9];
      char c2[5];
    } in_data[4][4];
    memset(&in_data[0][0], padding_byte_value, sizeof(in_data));
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        for (size_t k = 0; k < 3; ++k)
          in_data[i][j].c1[k] = (char)((i * 4 + j) * 3 + k);
        for (size_t k = 0; k < 2; ++k)
          in_data[i][j].i[k] = (int)((i * 4 + j) * 2 + k);
        for (size_t k = 0; k < 3; ++k)
          in_data[i][j].i16[k] = (int16_t)((i * 4 + j) * 3 + k);
        for (size_t k = 0; k < 9; ++k)
          for (size_t l = 0; l < 9; ++l)
            in_data[i][j].d[k][j] = (double)(((i * 4  + j) * 9 + k) * 9 + l);
        for (size_t k = 0; k < 5; ++k)
          in_data[i][j].c2[k] = (char)((i * 4 + j) * 5 + k);
      }
    }

    MPI_Datatype struct_mpi_ddt;
    {// STRUCT(3,{1,1,1},{...},{i_mpi_ddt,d_mpi_ddt,c_mpi_ddt})

      MPI_Datatype i_mpi_ddt, i16_mpi_ddt, d_mpi_ddt, c_mpi_ddt;
      { // {CONT(minimal_cont_count,DUP(MPI_INT)),
        //  CONT(3,MPI_INT16_T),
        //  VECTOR(3,3,9,MPI_DOUBLE),
        //  CONT(minimal_cont_count,MPI_CHAR)}

        { // CONT(minimal_cont_count,DUP(MPI_INT)
          MPI_Datatype dup_mpi_int;
          MPI_Type_dup(MPI_INT, &dup_mpi_int);
          enum {COUNT = minimal_cont_count};
          xt_mpi_call(MPI_Type_contiguous(COUNT, dup_mpi_int, &i_mpi_ddt),
                      MPI_COMM_WORLD);
          MPI_Type_free(&dup_mpi_int);
        }

        { // CONT(3,MPI_INT16_T)
          enum {COUNT = 3};
          MPI_Type_contiguous(COUNT, MPI_INT16_T, &i16_mpi_ddt);
        }

        { // VECTOR(3,3,9,MPI_DOUBLE)
          enum {COUNT = 3, BLOCKLENGTH = 3, STRIDE = 9};
          MPI_Type_vector(COUNT, BLOCKLENGTH, STRIDE, MPI_DOUBLE, &d_mpi_ddt);
        }

        { // CONT(0,MPI_CHAR)
          enum {COUNT = minimal_cont_count};
          MPI_Type_contiguous(COUNT, MPI_CHAR, &c_mpi_ddt);
        }
      }

      enum {COUNT = 4};
      static const int array_of_blocklengths[COUNT] = {1, 1, 1, 1};
      MPI_Aint base_address, i_address, i16_address, d_address, c_address;
      MPI_Get_address(in_data, &base_address);
      MPI_Get_address(in_data[0][0].i, &i_address);
      MPI_Get_address(in_data[0][0].i16, &i16_address);
      MPI_Get_address(&(in_data[0][0].d[3][3]), &d_address);
      MPI_Get_address(in_data[0][0].c2, &c_address);
      MPI_Aint array_of_displacements[COUNT] =
        { i_address - base_address,
          i16_address - base_address,
          d_address - base_address,
          c_address - base_address };
      const MPI_Datatype array_of_types[COUNT] =
        {i_mpi_ddt, i16_mpi_ddt, d_mpi_ddt, c_mpi_ddt};
      MPI_Type_create_struct(
        COUNT, (int *)array_of_blocklengths, array_of_displacements,
        (MPI_Datatype *)array_of_types, &struct_mpi_ddt);

      MPI_Type_free(&i_mpi_ddt);
      MPI_Type_free(&i16_mpi_ddt);
      MPI_Type_free(&d_mpi_ddt);
      MPI_Type_free(&c_mpi_ddt);
    }

    MPI_Datatype mpi_ddt;
    enum {COUNT = 4};
    static const int array_of_blocklengths[COUNT] = {1,2,3,4};
    static const int array_of_displacements[COUNT] = {3,7,9,12};
    MPI_Type_indexed(
      COUNT, (int *)array_of_blocklengths, (int *)array_of_displacements,
      struct_mpi_ddt, &mpi_ddt);
    MPI_Type_free(&struct_mpi_ddt);

    // force generation of ddt data
    if (xt_ddt_get_pack_size(mpi_ddt) == 0) PUT_ERR("");
    MPI_Datatype mpi_ddt_cpy;
    MPI_Type_dup(mpi_ddt, &mpi_ddt_cpy);

    size_t ref_pack_size
      = (1 + 2 + 3 + 4)
      * (minimal_cont_count * sizeof (int)
         + 3 * sizeof(int16_t)
         + 3 * 3 * sizeof(double)
         + minimal_cont_count);
    if (check_xt_ddt(
          mpi_ddt, in_data, sizeof(in_data), ref_pack_size))
      PUT_ERR("ERROR INDEXED(...)");

    if (check_xt_ddt(
          mpi_ddt_cpy, in_data, sizeof(in_data), ref_pack_size))
      PUT_ERR("ERROR DUP(INDEXED(...))");
  }

  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static int check_xt_ddt_off(
  MPI_Datatype mpi_ddt, const void *in_data_, int in_data_offset,
  size_t in_data_size, size_t ref_pack_size) {

  if (in_data_size > INT_MAX) {
    PUT_ERR("in_data_size > INT_MAX");
    return 0;
  }

  const unsigned char *in_data = in_data_;

  // check type of the MPI Datatype
  int num_ints, num_adds, num_dtypes, combiner;
  MPI_Type_get_envelope(
    mpi_ddt, &num_ints, &num_adds, &num_dtypes, &combiner);

  // if mpi_ddt is a derived MPI datatype
  if (combiner != MPI_COMBINER_NAMED)
    MPI_Type_commit(&mpi_ddt);

  if (xt_ddt_get_pack_size(mpi_ddt) != ref_pack_size) return 1;

  unsigned char *pack_data = xmalloc(ref_pack_size);
  unsigned char *out_data = xmalloc(in_data_size);
  unsigned char *ref_out_data = xmalloc(in_data_size);
  memset(ref_out_data, padding_byte_value, in_data_size);

  // generate reference data
  {
    /* can't use ref_pack_size for MPI buffer because implementations
     * are free to do whatever */
    int pack_size_for_mpi;
    MPI_Pack_size(1, mpi_ddt, MPI_COMM_WORLD, &pack_size_for_mpi);
    void *pack_buf = xmalloc((size_t)pack_size_for_mpi);
    int position = 0;
    MPI_Pack(
      CAST_MPI_SEND_BUF(in_data + in_data_offset),
      1, mpi_ddt, pack_buf, pack_size_for_mpi, &position, MPI_COMM_WORLD);
    position = 0;
    MPI_Unpack(
      pack_buf, pack_size_for_mpi, &position,
      ref_out_data + in_data_offset, 1, mpi_ddt, MPI_COMM_WORLD);
    free(pack_buf);
  }

  int err_flag = 0;

PragmaACC(data copyin(in_data[:in_data_size]) \
          create(pack_data[:ref_pack_size], out_data[:in_data_size]))
  {
#ifdef _OPENACC
    for (int data_on_gpu = 0; data_on_gpu < 2; ++data_on_gpu) {
      for (int buffer_on_gpu = 0; buffer_on_gpu < 2; ++buffer_on_gpu) {
#endif

        // initialise out data (also on device, if necessary)
        memset(out_data, padding_byte_value, in_data_size);
PragmaACC(update device(out_data[:in_data_size]) if (data_on_gpu))

        // use device pointers of in_data and out_data, if they
        // are on the device
#if defined _OPENACC && _OPENACC >= 201711
PragmaACC(host_data use_device(in_data, out_data) if(data_on_gpu))
#endif
        {
#if defined _OPENACC && _OPENACC < 201711
          void *in_tmp = in_data, *out_tmp = out_data;
          unsigned char *in_data = data_on_gpu ? acc_deviceptr(in_tmp) : in_tmp,
            *out_data = data_on_gpu ? acc_deviceptr(out_tmp) : out_tmp;
          void *tmp = pack_data;
          void *pack_data = buffer_on_gpu ? acc_deviceptr(tmp) : tmp;
#endif
          // use device pointer of pack_data, if it is on the device
#if defined _OPENACC && _OPENACC >= 201711
PragmaACC(host_data use_device(pack_data) if (buffer_on_gpu))
#endif
          {
            // pack in_data
            xt_ddt_pack(
              mpi_ddt, in_data + in_data_offset,
              pack_data);

            // unpack data
            xt_ddt_unpack(
              mpi_ddt, pack_data, out_data + in_data_offset);
          }
        }

        // copy data from device to host, if necessary
PragmaACC(update self(out_data[:in_data_size]) if (data_on_gpu))

        // compare out buffers
        err_flag |= memcmp(out_data, ref_out_data, in_data_size);

#ifdef _OPENACC
      } // buffer_on_gpu
    } // data_on_gpu
#endif
  }

  // clean up
  free(ref_out_data);
  free(out_data);
  free(pack_data);

  // if mpi_ddt is a derived MPI datatype
  if (combiner != MPI_COMBINER_NAMED) MPI_Type_free(&mpi_ddt);

  return err_flag;
}

static int check_xt_ddt(
  MPI_Datatype mpi_ddt, void const * in_data, size_t in_data_size,
  size_t ref_pack_size) {

  return check_xt_ddt_off(mpi_ddt, in_data, 0, in_data_size, ref_pack_size);
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
