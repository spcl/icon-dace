#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "yac.h"

int main(int argc, char** argv){
  if(argc != 2){
    fprintf(stderr, "Please provide the variable name as command-line argument!\n");
    return 1;
  }
  const char* field_name = argv[1];

  yac_cinit();

  int comp_id;
  yac_cdef_comp("yac_input_external", &comp_id);

  double* x_vertices = malloc(360*sizeof(*x_vertices));
  for(int i = 0; i<360;++i){
    x_vertices[i] = (M_PI/180.)*i;
  }
  double* y_vertices = malloc(181*sizeof(*y_vertices));
  for(int i = -90; i<=90;++i){
    y_vertices[i+90] = (M_PI/180.)*i;
  }
  int grid_id;
  int dims[] = {360, 181};
  int cyclic[] = {1,0};
  const char * grid_name = "yac_input_external_grid";
  yac_cdef_grid_reg2d ( grid_name, dims, cyclic,
                        x_vertices, y_vertices, &grid_id);

  int point_id;
  yac_cdef_points_reg2d(grid_id, dims, YAC_LOCATION_CORNER,
                        x_vertices, y_vertices, &point_id);
  free(x_vertices);
  free(y_vertices);

  yac_csync_def();
  const char* dt = yac_cget_field_timestep ( "yac_input_plugin",
                                             "comin_yac_plugin_grid", field_name);

  int field_id;
  yac_cdef_field(field_name, comp_id, &point_id, 1, 1,
                 dt, YAC_TIME_UNIT_ISO_FORMAT, &field_id);

  int interp_id;
  yac_cget_interp_stack_config(&interp_id);
  yac_cadd_interp_stack_config_nnn(interp_id, YAC_NNN_AVG, 1, 1.0);

  yac_cdef_couple(
    "yac_input_external", grid_name, field_name,
    "yac_input_plugin", "comin_yac_plugin_grid", field_name,
    dt, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
    interp_id, 1, 0);
  yac_cfree_interp_stack_config(interp_id);

  yac_cenddef();

  double* data = malloc(dims[0]*dims[1]*sizeof(*data));
  int info = YAC_ACTION_NONE;
  int ierror = 0;
  int t = 0;
  while(info != YAC_ACTION_OUT_OF_BOUND && info != YAC_ACTION_PUT_FOR_RESTART){
    // fill with artificial data
    for(int i = 0; i<dims[0]; ++i){
      for(int j = 0; j<dims[1]; j++){
        double lon = (M_PI/180.)*i;
        double lat = (M_PI/180.)*j;
        data[i*dims[1]+j] = cos(lat)*sin(lon + t*0.1);
      }
    }

    const char* time = yac_cget_field_datetime(field_id);
    yac_cput_ ( field_id, 1,data, &info, &ierror );
    fprintf(stderr, "datetime: %s: yac info: %d\n", time, info);
    assert(ierror == 0);
    ++t;
  }

  yac_cfinalize();
  return 0;
}
