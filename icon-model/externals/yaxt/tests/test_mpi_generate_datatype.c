/**
 * @file test_mpi_generate_datatype.c
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "ctest_common.h"

#define FIXED_MPI_ON_BLIZZARD

#ifdef FIXED_MPI_ON_BLIZZARD

static int
test_datatype_int(MPI_Datatype datatype, int recv_count, const int *send_data,
                  const int *ref_recv_data) {

  int recv_data[recv_count];
  MPI_Status    status;


  //datatype
  xt_mpi_call(MPI_Sendrecv(CAST_MPI_SEND_BUF(send_data), 1, datatype, 0, 0,
                           recv_data, recv_count, MPI_INT, 0, 0,
                           MPI_COMM_WORLD, &status),
              MPI_COMM_WORLD);

  bool mismatch = false;
  for (int i = 0; i < recv_count; ++i)
    mismatch |= (recv_data[i] != ref_recv_data[i]);

  return mismatch;
}

#else

static int
test_datatype_int(MPI_Datatype datatype, int recv_count, const int *send_data,
                  const int *ref_recv_data) {

  int recv_data[recv_count];
  MPI_Status    status;
  MPI_Request request;

  xt_mpi_call(MPI_Irecv(recv_data, recv_count, MPI_INT, 0, 0, MPI_COMM_WORLD,
                        &request), MPI_COMM_WORLD);

  xt_mpi_call(MPI_Send(CAST_MPI_SEND_BUF(send_data), 1, datatype, 0, 0,
                       MPI_COMM_WORLD), MPI_COMM_WORLD);

  xt_mpi_call(MPI_Wait(&request, &status), MPI_COMM_WORLD);

  bool mismatch = false;
  for (int i = 0; i < recv_count; ++i)
    mismatch |= (recv_data[i] != ref_recv_data[i]);

  return mismatch;
}

#endif

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // tests with empty data
    enum { count = 0 };
    static const int displacements[1] = { 0 };
    static const int blocklengths[1] = { 0 };
    MPI_Datatype old_type = MPI_INT;

    MPI_Datatype datatype
      = xt_mpi_generate_datatype(displacements, count, old_type,
                                 MPI_COMM_WORLD);

    if (datatype != MPI_DATATYPE_NULL)
      PUT_ERR("error in xt_mpi_generate_datatype (count == 0)\n");

    datatype = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                              count, old_type, MPI_COMM_WORLD);

    if (datatype != MPI_DATATYPE_NULL)
      PUT_ERR("error in xt_mpi_generatxt_mpi_generate_datatype_blocke_datatype"
              " (count == 0)\n");
  }

  { // tests with block_count == 1

    { // blocklength == 1 and displacement = 0
      enum { count = 1 };
      static const int displacements[count] = {0};
      static const int blocklengths[count] = {1};
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[3] = { 123, 234, 345 };
      static const int ref_recv_data[1] = { 123 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(displacements, count, old_type,
                                   MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype
        = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                         count, old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by"
                " xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // blocklength == 1 and displacement = 1
      enum { count = 1 };
      static const int displacements[count] = { 1 };
      static const int blocklengths[count] = { 1 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[3] = { 123, 234, 345 };
      static const int ref_recv_data[count] = { 234 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(displacements, count, old_type,
                                   MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype
        = xt_mpi_generate_datatype_block(displacements, blocklengths,
                                         count, old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);

    }

    { // blocklength == 2 and displacement = 0
      MPI_Datatype datatype;
      enum {
        element_count = 2,
        block_count = 1,
      };
      static const int element_displacements[element_count] = {0,1};
      static const int block_displacements[block_count] = {0};
      static const int blocklengths[block_count] = {2};
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[3] = { 123, 234, 345 };
      static const int ref_recv_data[element_count] = { 123, 234 };

      datatype = xt_mpi_generate_datatype(element_displacements, element_count,
                                          old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // blocklength == 2 and displacement = 1
      enum {
        element_count = 2,
        block_count = 1,
      };
      static const int element_displacements[element_count] = {1,2};
      static const int block_displacements[block_count] = {1};
      static int blocklengths[block_count] = {2};
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[3] = { 123, 234, 345 };
      static const int ref_recv_data[element_count] = { 234, 345 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths,
                                                block_count, old_type,
                                                MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }
  }

  { // tests with block_count == 3
    { // block length == {1,1,1} and block displacement = {0,2,4}
      enum {
        element_count = 3,
        block_count = 3,
      };
      static const int element_displacements[element_count] = { 0, 2, 4 };
      static const int block_displacements[block_count] = { 0, 2, 4 };
      static const int blocklengths[block_count] = { 1, 1, 1 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[5] = { 123, 234, 345, 456, 567 };
      static const int ref_recv_data[3] = { 123, 345, 567 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {0,3,6}
      enum {
        element_count = 6,
        block_count = 3,
      };
      static const int element_displacements[element_count]
        = { 0, 1, 3, 4, 6, 7 };
      static const int block_displacements[block_count] = { 0, 3, 6 };
      static const int blocklengths[block_count] = { 2, 2, 2 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
      static const int ref_recv_data[element_count] = { 0, 1, 3, 4, 6, 7 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {1,4,7}
      enum {
        element_count = 6,
        block_count = 3,
      };
      static const int element_displacements[element_count]
        = { 1, 2, 4, 5, 7, 8};
      static const int block_displacements[block_count] = { 1, 4, 7};
      static const int blocklengths[block_count] = { 2, 2, 2 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
      static const int ref_recv_data[element_count] = { 1, 2, 4, 5, 7, 8 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype
        = xt_mpi_generate_datatype_block(block_displacements,
                                         blocklengths, block_count,
                                         old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {7,4,1}
      enum {
        element_count = 6,
        block_count = 3,
      };
      static const int element_displacements[element_count]
        = { 7, 8, 4, 5, 1, 2 };
      static const int block_displacements[block_count] = { 7, 4, 1 };
      static const int blocklengths[block_count] = { 2, 2, 2 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
      static const int ref_recv_data[element_count] = { 7, 8, 4, 5, 1, 2 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,2,2} and block displacement = {0,3,7}
      enum {
        element_count = 6,
        block_count = 3,
      };
      static const int element_displacements[element_count]
        = { 0, 1, 3, 4, 7, 8 };
      static const int block_displacements[block_count] = { 0, 3, 7 };
      static const int blocklengths[block_count] = { 2, 2, 2 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
      static const int ref_recv_data[element_count] = { 0, 1, 3, 4, 7, 8 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }

    { // block length == {2,3,4} and block displacement = {0,3,7}
      enum {
        element_count = 9,
        block_count = 3,
      };
      static const int element_displacements[element_count]
        = { 0, 1, 3, 4, 5, 7, 8, 9, 10 };
      static const int block_displacements[block_count] = { 0, 3, 7 };
      static const int blocklengths[block_count] = { 2, 3, 4 };
      MPI_Datatype old_type = MPI_INT;

      static const int send_data[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
      static const int ref_recv_data[element_count]
        = { 0, 1, 3, 4, 5, 7, 8, 9, 10 };

      MPI_Datatype datatype
        = xt_mpi_generate_datatype(element_displacements, element_count,
                                   old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

      MPI_Type_free(&datatype);

      datatype = xt_mpi_generate_datatype_block(block_displacements,
                                                blocklengths, block_count,
                                                old_type, MPI_COMM_WORLD);

      if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
        PUT_ERR("error in datatype generated by "
                "xt_mpi_generate_datatype_block\n");

      MPI_Type_free(&datatype);
    }
  }


  { // check two encoded trivial vectors
    enum { element_count = 4 };
    static const int element_displacements[element_count] = { 2, 4, 7, 9 };
    MPI_Datatype old_type = MPI_INT;

    static const int send_data[]
      = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
          10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
          20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
    static const int ref_recv_data[element_count] = { 2, 4, 7, 9 };
    MPI_Datatype datatype
      = xt_mpi_generate_datatype(element_displacements, element_count,
                                 old_type, MPI_COMM_WORLD);

    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);

  }

  { // check two encoded non-trivial vectors
    enum { element_count = 8 };
    static const int element_displacements[element_count]
      = { 1, 4, 7, 10, 11, 14, 17, 20};
    MPI_Datatype old_type = MPI_INT;

    static const int send_data[]
      = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
          10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
          20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
    static const int ref_recv_data[element_count]
      = { 1, 4, 7, 10, 11, 14, 17, 20};
    MPI_Datatype datatype
      = xt_mpi_generate_datatype(element_displacements, element_count,
                                 old_type, MPI_COMM_WORLD);
    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);
  }

  { // check reverse block
    enum { element_count = 14 };
    static const int element_displacements[element_count]
      = {14,13,12,11,10,9,8,7,6,5,4,3,2,1};

    MPI_Datatype old_type = MPI_INT;

    static const int send_data[]
      = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
          10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
          20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
    static const int ref_recv_data[element_count]
      = { 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
    MPI_Datatype datatype
      = xt_mpi_generate_datatype(element_displacements, element_count,
                                 old_type, MPI_COMM_WORLD);
    if (test_datatype_int(datatype, element_count, send_data, ref_recv_data))
      PUT_ERR("error in datatype generated by xt_mpi_generate_datatype\n");

    MPI_Type_free(&datatype);

  }

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
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
