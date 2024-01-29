/**
   @file comin.h
   @brief C interface for the ICOM Community Interface

   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#ifndef COMIN_H
#define COMIN_H

#include <stdint.h>
#include <stdbool.h>
#include "comin_version.inc"
#include "comin_global.inc"

#ifdef __cplusplus
typedef bool _Bool;
extern "C" {
#endif

/// @defgroup c_interface C Interface
/// @{

  struct t_comin_var_descriptor{
    char name[MAX_LEN_VAR_NAME+1];
    int id;
  };

  enum ENTRY_POINT{
    EP_SECONDARY_CONSTRUCTOR = 1,
    EP_ATM_YAC_DEFCOMP_BEFORE,
    EP_ATM_YAC_DEFCOMP_AFTER,
    EP_ATM_YAC_SYNCDEF_BEFORE,
    EP_ATM_YAC_SYNCDEF_AFTER,
    EP_ATM_YAC_ENDDEF_BEFORE,
    EP_ATM_YAC_ENDDEF_AFTER,
    EP_ATM_INIT_FINALIZE,
    EP_ATM_TIMELOOP_BEFORE,
    EP_ATM_TIMELOOP_START,
    EP_ATM_TIMELOOP_END,
    EP_ATM_TIMELOOP_AFTER,
    EP_ATM_INTEGRATE_BEFORE,
    EP_ATM_INTEGRATE_START,
    EP_ATM_INTEGRATE_END,
    EP_ATM_INTEGRATE_AFTER,
    EP_ATM_WRITE_OUTPUT_BEFORE,
    EP_ATM_WRITE_OUTPUT_AFTER,
    EP_ATM_CHECKPOINT_BEFORE,
    EP_ATM_CHECKPOINT_AFTER,
    EP_ATM_ADVECTION_BEFORE,
    EP_ATM_ADVECTION_AFTER,
    EP_ATM_PHYSICS_BEFORE,
    EP_ATM_PHYSICS_AFTER,
    EP_ATM_NUDGING_BEFORE,
    EP_ATM_NUDGING_AFTER,
    EP_ATM_SURFACE_BEFORE,
    EP_ATM_SURFACE_AFTER,
    EP_ATM_TURBULENCE_BEFORE,
    EP_ATM_TURBULENCE_AFTER,
    EP_ATM_MICROPHYSICS_BEFORE,
    EP_ATM_MICROPHYSICS_AFTER,
    EP_ATM_CONVECTION_BEFORE,
    EP_ATM_CONVECTION_AFTER,
    EP_ATM_RADIATION_BEFORE,
    EP_ATM_RADIATION_AFTER,
    EP_ATM_RADHEAT_BEFORE,
    EP_ATM_RADHEAT_AFTER,
    EP_ATM_GWDRAG_BEFORE,
    EP_ATM_GWDRAG_AFTER,
    EP_FINISH,
    EP_DESTRUCTOR
  };

  enum ZAXIS {
    COMIN_ZAXIS_NONE         = 0,
    COMIN_ZAXIS_2D           = 1,
    COMIN_ZAXIS_3D           = 2,
    COMIN_ZAXIS_UNDEF        = 3,
  };

  enum VARACCESS_FLAG {
    FLAG_NONE         = 0,
    COMIN_FLAG_READ         = 1 << 1,
    COMIN_FLAG_WRITE        = 1 << 2,
    // note: the COMIN_FLAG_SYNCHRONIZED is unavailable (not yet implemented):
    // COMIN_FLAG_SYNCHRONIZED = 1 << 3
  };

  enum ERROR_CODE {
    COMIN_SUCCESS = 0,
    COMIN_INFO,
    COMIN_WARNING,
    COMIN_ERROR,
    COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_CALLBACK_COMPLETE,
    COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN,
    COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,
    COMIN_ERROR_DESCRDATA_FINALIZE,
    COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_METADATA_KEY_NOT_FOUND,
    COMIN_ERROR_METADATA_GET_INSIdE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_SETUP_FINALIZE,
    COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED,
    COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,
    COMIN_ERROR_PLUGIN_INIT_PRECISION,
    COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,
    COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED,
    COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,
    COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,
    COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,
    COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,
    COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,
    COMIN_ERROR_VARIABLE_NOT_TRACER_IN_TURBULENT_TRANSPORT,
    COMIN_ERROR_VARIABLE_NOT_TRACER_IN_CONVECTIVE_TRANSPORT,
    COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,
    COMIN_ERROR_FIELD_NOT_ALLOCATED,
    COMIN_ERROR_POINTER_NOT_ASSOCIATED,
    COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,
    COMIN_ERROR_FATAL
  };

  const int COMIN_DOMAIN_OUTSIDE_LOOP = -1;

  int     comin_setup_get_verbosity_level();

  int     comin_current_get_ep();
  int     comin_current_get_domain_id();
  int     comin_current_get_plugin_id();
  void    comin_current_get_plugin_name(char const ** val, int* len, int* ierr);
  void    comin_current_get_plugin_options(char const ** val, int* len, int* ierr);
  void    comin_current_get_plugin_comm(char const ** val, int* len, int* ierr);
  void    comin_current_get_datetime(char const ** val, int* len, int* ierr);

  int     comin_parallel_get_plugin_mpi_comm();
  int     comin_parallel_get_host_mpi_comm();
  int     comin_parallel_get_host_mpi_rank();

  void    comin_plugin_finish(const char* routine, const char* text);

  void    comin_var_request_add(struct t_comin_var_descriptor var_descriptor, _Bool lmodexclusive, int* ierr);
  void*   comin_var_get(int context_len, int* context, struct t_comin_var_descriptor var_descriptor,int flag);
  double* comin_var_get_ptr(void* handle);
  void    comin_var_get_shape(void* handle, int shape[5], int* ierr);
  void    comin_var_get_pos(void* handle, int* pos_jc, int* pos_jk, int* pos_jb, int* pos_jn, int* ierr);
  void    comin_var_get_ncontained(void* handle, int* ncontained, int* ierr);
  void*   comin_var_get_descr_list_head();
  void*   comin_var_get_descr_list_next(void* current);
  void    comin_var_get_descr_list_var_desc(void* current, struct t_comin_var_descriptor* var_desc, int* ierr);

  typedef void (*CALLBACK_PTR)();
  void    comin_callback_register(int entry_point_id, CALLBACK_PTR fct_ptr, int* ierr);
  void    comin_callback_get_ep_name(int iep, char out_ep_name[MAX_LEN_EP_NAME+1], int* ierr);

  int     comin_metadata_get_typeid(const char* key);
  void    comin_metadata_set_integer(struct t_comin_var_descriptor var_descriptor, const char* key, int val, int* ierr);
  void    comin_metadata_set_logical(struct t_comin_var_descriptor var_descriptor, const char* key, _Bool val, int* ierr);
  void    comin_metadata_set_real(struct t_comin_var_descriptor var_descriptor, const char* key, double val, int* ierr);
  void    comin_metadata_set_character(struct t_comin_var_descriptor var_descriptor, const char* key, char const * val, int* ierr);
  void    comin_metadata_get_integer(struct t_comin_var_descriptor var_descriptor, const char* key, int* val, int* ierr);
  void    comin_metadata_get_logical(struct t_comin_var_descriptor var_descriptor, const char* key, _Bool* val, int* ierr);
  void    comin_metadata_get_real(struct t_comin_var_descriptor var_descriptor, const char* key, double* val, int* ierr);
  void    comin_metadata_get_character(struct t_comin_var_descriptor var_descriptor, const char* key, char const ** val, int* len, int* ierr);

  double  comin_descrdata_get_timesteplength(int jg);
  int     comin_descrdata_get_index(int j);
  int     comin_descrdata_get_block(int j);
  void    comin_descrdata_get_cell_indices(int jg, int i_blk, int i_startblk, int i_endblk, int* i_startidx, int* i_endidx, int irl_start, int irl_end);
  int     comin_descrdata_get_cell_npromz(int jg);
  int     comin_descrdata_get_edge_npromz(int jg);
  int     comin_descrdata_get_vert_npromz(int jg);
  int     comin_descrdata_index_lookup_glb2loc_cell(int jg, int global_idx);
  void    comin_descrdata_get_simulation_interval_exp_start(char const ** val, int* len, int* ierr);
  void    comin_descrdata_get_simulation_interval_exp_stop(char const ** val, int* len, int* ierr);
  void    comin_descrdata_get_simulation_interval_run_start(char const ** val, int* len, int* ierr);
  void    comin_descrdata_get_simulation_interval_run_stop(char const ** val, int* len, int* ierr);



  /// returns version info.
  static inline void comin_setup_get_version(unsigned int* major, unsigned int* minor, unsigned int* patch) {
    (*major) = COMIN_VERSION_MAJOR;
    (*minor) = COMIN_VERSION_MINOR;
    (*patch) = COMIN_VERSION_PATCH;
  }

  /// Convenience operation for accessing 2D/3D fields.
  static inline double* comin_var_to_3d(void* handle){
    int ierr;
    int pos_jc, pos_jb, pos_jk, pos_jn;
    comin_var_get_pos(handle, &pos_jc, &pos_jb, &pos_jk, &pos_jn, &ierr);
    if(ierr != 0) comin_plugin_finish("comin_var_to_3d", "comin_var_get_pos failed");
    if(pos_jc != 0 || pos_jb != 1 || pos_jk != 2)
      comin_plugin_finish("comin_var_to_3d", "comin_var_to_3d only works for pos_jc, pos_jb, pos_jk = 1,2,3 (fortran dimension indices)");
    double* ptr = comin_var_get_ptr(handle);
    return ptr;
  }

  // ---------------------------------------------------------

  // The following **internal** struct is used to describe the
  // descrdata structure dynamically, mainly for use in the
  // python_adapter.
  struct comin_descrdata_property_t{
    const char* name;
    void* get_function;
    const char* datatype; // c datatype
    int ndims;
    bool has_jg;
    const struct comin_descrdata_property_t* subtypes;
  };

#ifdef __cplusplus
} // extern C
#endif

/* Header extension for get of grid data and domain routines generated by python script (comin_header_c_ext_descrdata_get_domain.h.py) in utils/.  */
#include "comin_header_c_ext_descrdata_query_domain.h"

/* Header extension for get of global data routines generated by python script (comin_header_c_ext_descrdata_get_global.h.py) in utils/.  */
#include "comin_header_c_ext_descrdata_query_global.h"


/// @}

#endif
