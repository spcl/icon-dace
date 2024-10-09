// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "dist_merge.h"
#include "utils_common.h"

static size_t get_pack_size(
  size_t count, unsigned char * array, size_t element_size, MPI_Comm comm,
  size_t(*element_get_pack_size)(void * element, MPI_Comm comm)) {

  if (count == 0) return 0;

  int count_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &count_pack_size), comm);

  size_t array_pack_size = 0;
  for (size_t i = 0; i < count; ++i)
    array_pack_size +=
      element_get_pack_size((void*)(array + i * element_size), comm);

  return (size_t)count_pack_size + array_pack_size;
}

static void pack(
  size_t count, unsigned char * array, size_t element_size, void * buffer,
  int buffer_size, int * position, MPI_Comm comm,
  void(*element_pack)(
    void * element, void * buffer, int buffer_size,
    int * position, MPI_Comm)) {

  int count_int = (int)count;
  yac_mpi_call(
    MPI_Pack(
      &count_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  for (size_t i = 0; i < count; ++i)
    element_pack(
      (void*)(array + i * element_size), buffer, buffer_size, position, comm);
}

static void unpack(
  void * buffer, int buffer_size, int * position, size_t * count,
  unsigned char ** array, size_t element_size, MPI_Comm comm,
  void(*element_unpack)(
    void * buffer, int buffer_size, int * position,
    void * element, MPI_Comm comm)) {

  int count_int;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &count_int, 1, MPI_INT, comm), comm);

  *count = (size_t)count_int;
  *array = xmalloc(*count * element_size);

  for (size_t i = 0; i < *count; ++i)
    element_unpack(
      buffer, buffer_size, position, (void*)(*array + i * element_size), comm);
}

void yac_dist_merge(
  size_t * count, void ** array, size_t element_size, MPI_Comm comm,
  struct yac_dist_merge_vtable * vtable, size_t ** idx_old_to_new) {

  int rank;
  MPI_Comm_rank(comm, &rank);

  // initialise
  unsigned char * input = *array;
  size_t input_len = *count;
  size_t* idx = NULL;
  if (idx_old_to_new) {
    *idx_old_to_new = xmalloc(input_len * sizeof(**idx_old_to_new));
    idx = xmalloc(input_len * sizeof(*idx));
    for(size_t i = 0;i<input_len;++i) idx[i] = i;
  }
  unsigned char * arr_new = NULL;
  size_t len_new = 0;

  YAC_ASSERT(
    input_len < INT_MAX, "ERROR(yac_dist_merge): too many elements");

  // sort input data so that the processing order on all processes is identical
  if (idx)
    yac_qsort_index(input, input_len, element_size, vtable->compare, idx);
  else
    qsort((void*)input, input_len, element_size, vtable->compare);

  // loop until no process has any unsychronised elements left
  void * buffer = NULL;
  while(1){

    // determine pack size of all remaining local data
    size_t pack_size =
      get_pack_size(
        input_len, input, element_size, comm, vtable->get_pack_size);
    YAC_ASSERT(
      pack_size <= LONG_MAX, "ERROR(yac_dist_merge): packing size too big");

    // determine rank with most amount of data
    struct {
      long pack_size;
      int rank;
    } data_pair = {.pack_size = (long)pack_size, .rank = rank};
    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &data_pair, 1, MPI_LONG_INT, MPI_MAXLOC, comm), comm);

    // if there is no more data to exchange, the sychronisation is finished
    if (data_pair.pack_size == 0) break;

    // allocate the buffer according to the processes with the most amount
    // of unsychronised data
    pack_size = (size_t)data_pair.pack_size;
    if (!buffer) buffer = xmalloc(pack_size);
    int position = 0;

    // the process with the most amount of data packs it
    if(data_pair.rank == rank)
      pack(
        input_len, input, element_size, buffer, pack_size, &position, comm,
        vtable->pack);

    // broadcast and unpack data (this only contains the basic information
    // required to identify an element)
    yac_mpi_call(
      MPI_Bcast(buffer, pack_size, MPI_PACKED, data_pair.rank, comm), comm);
    unsigned char * recved = NULL;
    size_t num_recved;
    position = 0;
    unpack(
      buffer, pack_size, &position, &num_recved, &recved, element_size, comm,
      vtable->unpack);

    // copy the received elements into the result array
    arr_new = xrealloc(arr_new, (len_new + num_recved)*element_size);
    memcpy(
      arr_new + len_new*element_size, (void*)recved, num_recved*element_size);
    free(recved);

    // merge received elements into the result array
    size_t input_idx, input_len_new, i = len_new;
    len_new += num_recved;
    for(input_idx = 0, input_len_new = 0; i < len_new; ++i) {
      void* recved_element = arr_new + i*element_size;

      // search for matching element in input list until an element is
      // found, which is "bigger" (as defined by the compare function) or the
      // end of the input list was reached
      int cmp = 0;
      void * input_element = NULL;
      while ((input_idx < input_len) &&
             (((cmp =
                  vtable->compare(
                    ((input_element = input + input_idx * element_size)),
                    recved_element))) < 0)) {

        // elements from the input list that matched with a received element
        // are removed from the input list
        if (input_idx != input_len_new) {
          memcpy(
            input + input_len_new * element_size,
            input_element, element_size);
          if (idx) idx[input_len_new] = idx[input_idx];
        }
        ++input_len_new;
        ++input_idx;
      }

      // if the end of the local input list was reached, no more merging is required
      if (input_idx == input_len) break;

      // if a matching element was found in the input list
      if (!cmp) {

        // merge input list element with received element
        if (vtable->merge)
          vtable->merge(recved_element, input_element, comm);
        // free the element from the input list
        vtable->free_data(input_element);
        // keep track on where the original element is in the result list
        if (idx_old_to_new) (*idx_old_to_new)[idx[input_idx]] = i;
        // upate input list idx
        input_idx++;

      // if no matching element was found in the input list
      } else if (vtable->merge) {
        // since the merge operation can potentially be collective, we have
        // to call it, even if no matching element was found
        vtable->merge(recved_element, NULL, comm);
      }
    }

    // if the end of the input list was reched:
    //   process remaining received elements
    if (vtable->merge)
      for(; i < len_new; ++i)
        vtable->merge(arr_new + i*element_size, NULL, comm);

    // if all received elements where processed:
    //   compress remaining elements in input list
    for (; input_idx < input_len; ++input_idx, ++input_len_new) {
      if (input_idx != input_len_new) {
        memcpy(
          input + input_len_new * element_size,
          input + input_idx * element_size, element_size);
        if (idx) idx[input_len_new] = idx[input_idx];
      }
    }
    // update length of input list
    input_len = input_len_new;
  }
  free(buffer);
  free(input);
  free(idx);
  *array = arr_new;
  *count = len_new;
}
