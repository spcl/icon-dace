/**
 * @file tr_gp2fs.c
 *
 * @copyright Copyright  (C)  2012 Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
 *         Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
 *             Moritz Hanke <hanke@dkrz.de>
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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#include <yaxt.h>

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))
//#define DEBUG
//#define CHECK_RESULTS

enum initializationType {
  INIT_NONE,
  INIT_SIMPLE,
  INIT_IDXLIST,
};

struct config {

  enum initializationType initialize;
   int x, y, z; // global data size (default 1024 * 512 * 1)

   struct {
      int x_len, y_len, z_len;        //local data size
      int local_size[3];
      Xt_int local_start[3];
   } src, dst; //!< source and target decomposion
};

enum direction { SOURCE = 0, DESTINATION = 1};

#ifdef DEBUG
static void
print_index(Xt_idxlist idxsection, int rank, enum direction dir) {
  static const char dirMsg[2][4] = { "SRC", "DST" };
  int num_idx = xt_idxlist_get_num_indices(idxsection);
  for (int p = 0; p < num_idx; p ++){
    Xt_int value;
    xt_idxlist_get_index_at_position (idxsection, p, &value);
    printf ("rank: %d , INDEX_%s: %"XT_INT_FMT" \n", rank, dirMsg[dir], value);
  }
}
#endif

static void
calc_blocks(int nprocs, int dim2, int dim1, int *dim2_block, int *dim1_block) {
  // temp num of blocks in x direction
  int tmp_dim1 = -1;
  // temp prop to compare last and actual iteration
  double tmp_prop = -1.0;

  int snprocs = (int)(sqrt(nprocs));

  // calculate global data proportion by division of bigger number by smaller number, so prop <= 1
  double prop = (dim2 < dim1)
    ? (double)dim2/(double)dim1 : (double)dim1/(double)dim2;

  //
  for (int f = snprocs; f > 0; f-- ){
    if (nprocs % f == 0){
      *dim1_block = nprocs / f;
      *dim2_block = f;

      double actu = fabs (prop - (double)*dim2_block/(double)*dim1_block);
      double past = fabs  (prop - tmp_prop) ;

      //compare if actual prop better then last. If, make it to actual and do next iteration, otherwise take last and break
      if (actu < past ){
        tmp_prop = (double)*dim2_block / (double)*dim1_block;
        tmp_dim1 = *dim1_block;
      } else {
        *dim1_block = tmp_dim1;
        *dim2_block = nprocs / tmp_dim1;
        break;
      }
    }
  }

  // adapt proportion to global data
  if(dim2 > dim1){
    *dim1_block = *dim2_block;
    *dim2_block = tmp_dim1;
  }
}

static void
init_data(int *data_array, int array_len, int j){
  // 0 = source
  if (j == 0)
    for (int i = 0; i < array_len; i++) {
      data_array[i] = 100 + i;
    }
  // 1 = dest
  else
    for (int i = 0; i < array_len; i++) {
      data_array[i] = -1;
    }
}

static void
init_data_with_index_list(int *data_array, int array_len,
                          Xt_idxlist idxlist) {

  if (xt_idxlist_get_num_indices(idxlist) != array_len) {

    fputs("error in init_data_with_indices: array length and index"
          " list length does not match\n", stderr);
    exit(EXIT_FAILURE);
  }

#ifdef CHECK_RESULTS
  for (int i = 0; i < array_len; ++i) {
    xt_idxlist_get_index_at_position(idxlist, i, data_array+i);
  }
#else
  memset(data_array, 0, sizeof(*data_array) * (size_t)array_len);
#endif
}

static int
check_data_against_index_list(int *data_array,
                              int array_len,
                              Xt_idxlist idxlist) {

#ifdef CHECK_RESULTS
  if (xt_idxlist_get_num_indices(idxlist) != array_len) {

    fputs("error in check_data_against_index_list: array lenght and"
          " index list length does not match\n", stderr);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < array_len; ++i) {

    Xt_int index;
    xt_idxlist_get_index_at_position(idxlist, i, &index);
    if (data_array[i] != index) {

      fputs("error in check_data_against_index_list: data does not match\n",
            stderr);
      return 1;
    }
  }
#else
  if (array_len < 0 || !data_array || !idxlist)
    return 1;
#endif
  return 0;
}

static void
calc_config(struct config *conf, int nprocs) {
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  // --- SRC
  // source num blocks in global structure
  int src_y_block = 1;
  int src_x_block = 1;

  // calculate num of blocks regarding global data proportion
  calc_blocks (nprocs, conf->y, conf->x, &src_y_block, &src_x_block);

  // source local data size in x, y and z direction - num of elements

  conf->src.z_len = conf->z;
  conf->src.y_len = conf->y / src_y_block + ((rank / src_x_block) < (conf->y % src_y_block));
  conf->src.x_len = conf->x / src_x_block + ((rank % src_x_block) < (conf->x % src_x_block));

  // --- DEST
  //destination num blocks in global structure
  int dst_z_block = 1;
  int dst_y_block = 1;
  //int dst_x_block = 1;

  calc_blocks (nprocs, conf->z, conf->y, &dst_z_block, &dst_y_block);

  //calculate actually length of local destination data array
  conf->dst.z_len = conf->z / dst_z_block + ((rank / dst_y_block) < (conf->z % dst_z_block));
  conf->dst.y_len = conf->y / dst_y_block + ((rank % dst_y_block) < (conf->y % dst_y_block));
  conf->dst.x_len = conf->x;

  // calculate parameters for source section
  conf->src.local_size[0] = conf->z;
  conf->src.local_size[1] = conf->src.y_len;
  conf->src.local_size[2] = conf->src.x_len;
  int src_block_pos_y = rank / src_x_block;
  int src_block_pos_x = rank % src_x_block;

  // how many blocks get an extra row or column (source_dim_blocks_bonus)
  int src_y_bb = conf->y % src_y_block;
  int src_x_bb = conf->x % src_x_block;

  conf->src.local_start[0] = 0;
  conf->src.local_start[1] = (Xt_int)((conf->y/ src_y_block) * src_block_pos_y
                                      + MIN(src_block_pos_y, src_y_bb));
  conf->src.local_start[2] = (Xt_int)((conf->x/ src_x_block) * src_block_pos_x
                                      + MIN(src_block_pos_x, src_x_bb));

  //dst
  conf->dst.local_size[0] = conf->dst.z_len;
  conf->dst.local_size[1] = conf->dst.y_len;
  conf->dst.local_size[2] = conf->x;

  int dst_block_pos_z = rank / dst_y_block;
  int dst_block_pos_y = rank % dst_y_block;

  // how many blocks get an extra row or column (source_dim_blocks_bonus)
  int dst_z_bb = conf->z % dst_z_block;
  int dst_y_bb = conf->y % dst_y_block;

  conf->dst.local_start[0] = (Xt_int)((conf->z/ dst_z_block) * dst_block_pos_z
                                      + MIN(dst_block_pos_z, dst_z_bb));
  conf->dst.local_start[1] = (Xt_int)((conf->y/ dst_y_block) * dst_block_pos_y
                                      + MIN(dst_block_pos_y, dst_y_bb));
  conf->dst.local_start[2] = 0;
}


static void
get_config(struct config *conf, int argc, char** argv, int nprocs){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int map = 0;
  int par = 0;

  static const struct {
    const char str[8];
    enum initializationType type;
  } initMap[] = {
    {  "none", INIT_NONE },
    {  "simple", INIT_SIMPLE },
    {  "idxlist", INIT_IDXLIST },
  };
  conf->initialize = INIT_NONE;
  enum { initMapSize = sizeof (initMap) / sizeof (initMap[0]) };
  while(( par = getopt(argc, argv, "x:y:z:i:")) != -1) {
    switch(par) {
      case 'x':
        conf->x = atoi(optarg);
        map |= 1;
        break;
      case 'y':
        conf->y = atoi(optarg);
        map |= 2;
        break;
      case 'z':
        conf->z = atoi(optarg);
        map |= 4;
        break;
      case '?':
        //printf("No arguments \n");
        break;
      case 'i':
        {
          size_t i;
          for (i = 0; i < initMapSize; ++i)
            if (!strcmp(optarg, initMap[i].str))
              break;
          if (i < initMapSize)
            conf->initialize = initMap[i].type;
          else {
            fprintf(stderr, "unknown type of initialization requested: %s\n",
                    optarg);
            fputs("available types: ", stderr);
            for (i = 0; i < initMapSize; ++i)
              fprintf(stderr, "%s%s", initMap[i].str,
                      i < initMapSize - 1 ? ", " : "\n");
          }
        }
        break;
    }
  }

  if (argc == 1) {
    if(rank == 0)
      printf("No arguments, set default: x=1024, y=512, z=1 \n");
    conf->x = 1024;
    conf->y = 512;
    conf->z = 1;
  } else if((map != 7) || (optind != 7 )) {
    if (rank == 0)
      printf(" Option fail \n Example: ./tr_gp2fs -x 10 -y 20 -z 100 \n Abort \n");
    exit(EXIT_FAILURE);
  }

  calc_config (conf, nprocs);
}

int main(int argc, char** argv) {

  //init mpi

  MPI_Init(NULL, NULL);

  xt_initialize (MPI_COMM_WORLD);

  int nprocs, rank;
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  struct config conf;
  int ret_val;

  {
    get_config(&conf, argc, argv, nprocs);

    // global param
    enum { num_dims = 3 };
    Xt_int global_size[3] = {(Xt_int)conf.z, (Xt_int)conf.y, (Xt_int)conf.x};

    // -- SOURCE index list by section
    Xt_int start = 0;

    //create source index section
    Xt_idxlist src_idxsection
      = xt_idxsection_new(start, num_dims, global_size,
                          conf.src.local_size, conf.src.local_start);

#ifdef DEBUG
    // 0 = soruce
    print_index(src_idxsection, rank, SOURCE);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // --- DEST-ination index list by section
    Xt_idxlist dst_idxsection
      = xt_idxsection_new(start, num_dims, global_size,
                          conf.dst.local_size, conf.dst.local_start);

    // 1 = dest
#ifdef DEBUG
    print_index(dst_idxsection, rank, DESTINATION);
#endif
    // xmap
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxsection, dst_idxsection, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist
      = xt_redist_p2p_new(xmap, MPI_INT);

    // --- local data array
    int *src_array, *dst_array;

    //allocate memory for source
    int src_len = conf.src.z_len * conf.src.y_len * conf.src.x_len;
    assert(src_len > 0);
    src_array = malloc(sizeof(*src_array)* (size_t)src_len);

    // allocate memory of actuall local data size
    int dst_len = conf.dst.z_len * conf.dst.y_len * conf.dst.x_len;
    assert(dst_len > 0);
    dst_array = calloc((size_t)dst_len, sizeof(dst_array));

    // init arrays with fake data
   switch (conf.initialize) {
   case INIT_NONE:
     break;
   case INIT_SIMPLE:
     init_data(src_array, src_len, 0);
     init_data(dst_array, dst_len, 1);
     break;
   case INIT_IDXLIST:
     init_data_with_index_list(src_array, src_len, src_idxsection);
     break;
   }

    //array poiter, especially necessary for data array number > 1
    int* src_array_p = src_array;
    int* dst_array_p = dst_array;

    //Exchange
    xt_redist_s_exchange1(redist, src_array_p, dst_array_p);

    ret_val = check_data_against_index_list(dst_array, dst_len, dst_idxsection);

    //clean up
    free (src_array);
    free (dst_array);

    xt_idxlist_delete(dst_idxsection);
    xt_idxlist_delete(src_idxsection);
    xt_xmap_delete(xmap);
    xt_redist_delete(redist);
  }

  MPI_Finalize();

  return ret_val;
}

/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
