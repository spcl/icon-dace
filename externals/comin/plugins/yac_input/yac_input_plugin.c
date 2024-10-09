#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <comin.h>
#include <mpi.h>
#include <yac.h>

int comp_id;
int field_id;
void* comin_var;
struct t_comin_var_descriptor descr;

void secondary_constructor();
void defcomp();
void create_grid_and_field();
void recv_field();

void comin_main(){
  const char* options = NULL;
  int options_len = 0;
  comin_current_get_plugin_options(&options, &options_len);
  strncpy(descr.name, options, options_len);
  descr.id = 1;


  comin_var_request_add(descr, false);

  comin_metadata_set_integer(descr, "zaxis_id", COMIN_ZAXIS_2D);

  comin_callback_register(EP_SECONDARY_CONSTRUCTOR, secondary_constructor);
  comin_callback_register(EP_ATM_YAC_DEFCOMP_BEFORE, defcomp);
  comin_callback_register(EP_ATM_YAC_SYNCDEF_BEFORE, create_grid_and_field);
  comin_callback_register(EP_ATM_TIMELOOP_START, recv_field);
}

void secondary_constructor(){
  int ep = EP_ATM_TIMELOOP_START;
  comin_var = comin_var_get(1, &ep, descr, COMIN_FLAG_WRITE);
  if(comin_var == NULL)
    comin_plugin_finish("yac_input_plugin", "Failed to get variable");
}

void defcomp(){
  int instance_id = comin_descrdata_get_global_yac_instance_id();
  yac_cpredef_comp_instance(instance_id, "yac_input_plugin", &comp_id);
  fprintf(stderr, "comp_id: %d\n", comp_id);
}

void create_grid_and_field(){
  int n_patch_cells = comin_descrdata_get_domain_cells_ncells(1);
  int n_patch_verts = comin_descrdata_get_domain_verts_nverts(1);

  int* num_vertices_per_cell = malloc(n_patch_cells*sizeof(*num_vertices_per_cell));
  for(int i=0; i<n_patch_cells; ++i)
    num_vertices_per_cell[i] = 3;

  int vlon_size[2];
  double* vlon;
  comin_descrdata_get_domain_verts_vlon(1, &vlon, vlon_size);
  int vlat_size[2];
  double* vlat;
  comin_descrdata_get_domain_verts_vlat(1, &vlat, vlat_size);

  // yac needs radiants
  double* vlon_rad = malloc(n_patch_verts*sizeof(*vlon_rad));
  double* vlat_rad = malloc(n_patch_verts*sizeof(*vlat_rad));
  for(int i = 0; i<n_patch_verts; ++i){
    vlon_rad[i] = vlon[i]/180.*M_PI;
    vlat_rad[i] = vlat[i]/180.*M_PI;
  }

  int vertex_idx_size[3]; // nproma, nblks, 3
  int* vertex_idx, *vertex_blk;
  comin_descrdata_get_domain_cells_vertex_idx(1, &vertex_idx, vertex_idx_size);
  comin_descrdata_get_domain_cells_vertex_blk(1, &vertex_blk, vertex_idx_size);
  int nproma = vertex_idx_size[0];

  int* vertex_idx_transposed = malloc(3*n_patch_cells*sizeof(*vertex_idx_transposed));
  for(int i=0;i<n_patch_cells; ++i)
    for(int k=0;k<3;++k){
      vertex_idx_transposed[k + 3*i] = vertex_idx[i + k*vertex_idx_size[0]*vertex_idx_size[1]] -1
        + nproma*(vertex_blk[i + k*vertex_idx_size[0]*vertex_idx_size[1]] -1);
    }

  int grid_id;
  yac_cdef_grid_unstruct ( "comin_yac_plugin_grid",
                           n_patch_verts,
                           n_patch_cells,
                           num_vertices_per_cell,
                           vlon_rad,
                           vlat_rad,
                           vertex_idx_transposed,
                           &grid_id);
  free(vlon_rad);
  free(vlat_rad);
  free(vertex_idx_transposed);

  int clon_size[2];
  double* clon;
  comin_descrdata_get_domain_cells_clon(1, &clon, clon_size);
  int clat_size[2];
  double* clat;
  comin_descrdata_get_domain_cells_clat(1, &clat, clat_size);

  // yac needs radiants
  double* clon_rad = malloc(n_patch_cells*sizeof(*clon_rad));
  double* clat_rad = malloc(n_patch_cells*sizeof(*clat_rad));
  for(int i = 0; i<n_patch_cells; ++i){
    clon_rad[i] = clon[i]/180.*M_PI;
    clat_rad[i] = clat[i]/180.*M_PI;
  }

  int point_id;
  yac_cdef_points_unstruct ( grid_id, n_patch_cells, YAC_LOCATION_CELL,
                             clon_rad, clat_rad, &point_id );

  free(clon_rad);
  free(clat_rad);

  double dt = comin_descrdata_get_timesteplength(1);

  char dt_str[16];
  sprintf(dt_str, "%d", ((int)dt*1000));
  yac_cdef_field ( descr.name, comp_id, &point_id, 1, 1,
                   dt_str, YAC_TIME_UNIT_MILLISECOND, &field_id);

  free(num_vertices_per_cell);
}

void recv_field(){
  int info, ierror;
  double* ptr = comin_var_get_ptr(comin_var);
  yac_cget_(field_id, 1, ptr, &info, &ierror);
}
