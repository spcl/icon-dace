// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DIST_MERGE_H
#define DIST_MERGE_H

#include "yac_mpi_common.h"

// The distribute merge is used by the coupling configuration in order to
// synchronise the configuration between all processes.

struct yac_dist_merge_vtable {

  /**
   * Determines pack size of an element.
   * The data to be packed only has to include the basic information that
   * is required to the element. Additional data is exchanged at a later stage.
   */
  size_t(*get_pack_size)(void * element, MPI_Comm comm);

  /**
   * Packs an element (only basic information).
   */
  void(*pack)(
    void * element, void * buffer, int buffer_size, int * position, MPI_Comm);

  /**
   * Unpacks an element (only basic information). The non-basic data of
   * the element is also to be initialised by this routine using dummy values.
   */
  void(*unpack)(
    void * buffer, int buffer_size, int * position, void * element,
    MPI_Comm comm);

  /**
   * Compares the basic information of two elements.
   */
  int(*compare)(void const * a, void const * b);

  /**
   * Merges two elements that have the same basic information.
   * The first argument contains an element that was initialised
   * by \c unpack and only contains the basic information. The second
   * argument contains a local element, which can contain additional
   * information. If a process does not have a matching element, local
   * will be \c NULL. The user is responsible for synchronising the
   * additional information between all processes in \c comm and add this
   * information to the element in the first argument. He also has to check
   * the consistency of the additional information between the processes.
   * The second argument will be freed after this call.
   *
   * This pointer can be zero, if an element only contains basic
   * information.
   */
  void(*merge)(void * to, void * from, MPI_Comm comm);

  /**
   * Frees a single element.
   */
  void(*free_data)(void * element);
};

/**
 * Distributes and merges an array across all processes in comm. After the call,
 * the array will be the same on all processes.
 * @param[inout] count          number of elements in the array
 * @param[inout] array          array of elements to by synchronised
 * @param[in]    element_size   size of a single element
 * @param[in]    comm           communicator
 * @param[in]    vtable         virtual method table that is to be used for
 *                              the synchronisation
 * @param[out]   idx_old_to_new if idx_old_to_new != NULL this routine will
 *                              allocate an array that has the original size
 *                              of the input array. After the call the i'th
 *                              element of the input array will be at position
 *                              idx_old_to_new[i] in the new array.
 * @remark This routine can be applied recursively: The \c merge function may
 *         itself call \c yac_dist_merge to synchronise some of the additional
 *         information.
 */
void yac_dist_merge(
  size_t * count, void ** array, size_t element_size, MPI_Comm comm,
  struct yac_dist_merge_vtable * vtable, size_t ** idx_old_to_new);

#endif // DIST_MERGE_H
