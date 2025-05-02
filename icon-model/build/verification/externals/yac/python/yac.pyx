# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

import cython
from enum import Enum
import numpy as _np
import logging
from types import coroutine

from libc.stdlib cimport malloc, free

cdef import from "<mpi.h>" nogil:
  ctypedef struct _mpi_comm_t
  ctypedef _mpi_comm_t* MPI_Comm
  int MPI_Comm_c2f(MPI_Comm)
  MPI_Comm MPI_Comm_f2c(int)

cdef extern from "Python.h":
    int Py_AtExit(void (*)())

cdef extern from "yac.h":
  cdef const int _LOCATION_CELL "YAC_LOCATION_CELL"
  cdef const int _LOCATION_CORNER "YAC_LOCATION_CORNER"
  cdef const int _LOCATION_EDGE "YAC_LOCATION_EDGE"

  cdef const int _EXCHANGE_TYPE_NONE "YAC_EXCHANGE_TYPE_NONE"
  cdef const int _EXCHANGE_TYPE_SOURCE "YAC_EXCHANGE_TYPE_SOURCE"
  cdef const int _EXCHANGE_TYPE_TARGET "YAC_EXCHANGE_TYPE_TARGET"

  cdef const int _ACTION_NONE "YAC_ACTION_NONE"
  cdef const int _ACTION_REDUCTION "YAC_ACTION_REDUCTION"
  cdef const int _ACTION_COUPLING "YAC_ACTION_COUPLING"
  cdef const int _ACTION_GET_FOR_RESTART "YAC_ACTION_GET_FOR_RESTART"
  cdef const int _ACTION_PUT_FOR_RESTART "YAC_ACTION_PUT_FOR_RESTART"
  cdef const int _ACTION_OUT_OF_BOUND "YAC_ACTION_OUT_OF_BOUND"

  cdef const int _REDUCTION_TIME_NONE "YAC_REDUCTION_TIME_NONE"
  cdef const int _REDUCTION_TIME_ACCUMULATE "YAC_REDUCTION_TIME_ACCUMULATE"
  cdef const int _REDUCTION_TIME_AVERAGE "YAC_REDUCTION_TIME_AVERAGE"
  cdef const int _REDUCTION_TIME_MINIMUM "YAC_REDUCTION_TIME_MINIMUM"
  cdef const int _REDUCTION_TIME_MAXIMUM "YAC_REDUCTION_TIME_MAXIMUM"

  cdef const int _CALENDAR_NOT_SET "YAC_CALENDAR_NOT_SET"
  cdef const int _PROLEPTIC_GREGORIAN "YAC_PROLEPTIC_GREGORIAN"
  cdef const int _YEAR_OF_365_DAYS "YAC_YEAR_OF_365_DAYS"
  cdef const int _YEAR_OF_360_DAYS "YAC_YEAR_OF_365_DAYS"

  cdef const int _TIME_UNIT_MILLISECOND "YAC_TIME_UNIT_MILLISECOND"
  cdef const int _TIME_UNIT_SECOND "YAC_TIME_UNIT_SECOND"
  cdef const int _TIME_UNIT_MINUTE "YAC_TIME_UNIT_MINUTE"
  cdef const int _TIME_UNIT_HOUR "YAC_TIME_UNIT_HOUR"
  cdef const int _TIME_UNIT_DAY "YAC_TIME_UNIT_DAY"
  cdef const int _TIME_UNIT_MONTH "YAC_TIME_UNIT_MONTH"
  cdef const int _TIME_UNIT_YEAR "YAC_TIME_UNIT_YEAR"
  cdef const int _TIME_UNIT_ISO_FORMAT "YAC_TIME_UNIT_ISO_FORMAT"

  cdef const int _AVG_ARITHMETIC "YAC_AVG_ARITHMETIC"
  cdef const int _AVG_DIST "YAC_AVG_DIST"
  cdef const int _AVG_BARY "YAC_AVG_BARY"

  cdef const int _NCC_AVG "YAC_NCC_AVG"
  cdef const int _NCC_DIST "YAC_NCC_DIST"

  cdef const int _NNN_AVG "YAC_NNN_AVG"
  cdef const int _NNN_DIST "YAC_NNN_DIST"
  cdef const int _NNN_GAUSS "YAC_NNN_GAUSS"
  cdef const int _NNN_RBF "YAC_NNN_RBF"

  cdef const int _CONSERV_DESTAREA "YAC_CONSERV_DESTAREA"
  cdef const int _CONSERV_FRACAREA "YAC_CONSERV_FRACAREA"

  cdef const int _SPMAP_AVG "YAC_SPMAP_AVG"
  cdef const int _SPMAP_DIST "YAC_SPMAP_DIST"

  cdef const int _SPMAP_NONE "YAC_SPMAP_NONE"
  cdef const int _SPMAP_SRCAREA "YAC_SPMAP_SRCAREA"
  cdef const int _SPMAP_INVTGTAREA "YAC_SPMAP_INVTGTAREA"
  cdef const int _SPMAP_FRACAREA "YAC_SPMAP_FRACAREA"

  cdef const int _CONFIG_OUTPUT_FORMAT_YAML "YAC_CONFIG_OUTPUT_FORMAT_YAML"
  cdef const int _CONFIG_OUTPUT_FORMAT_JSON "YAC_CONFIG_OUTPUT_FORMAT_JSON"

  cdef const int _CONFIG_OUTPUT_SYNC_LOC_DEF_COMP "YAC_CONFIG_OUTPUT_SYNC_LOC_DEF_COMP"
  cdef const int _CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF "YAC_CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF"
  cdef const int _CONFIG_OUTPUT_SYNC_LOC_ENDDEF "YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF"

  void yac_cinit ()
  void yac_cinit_instance ( int * yac_instance_id )
  void yac_cinit_comm (MPI_Comm comm )
  void yac_cinit_comm_instance (MPI_Comm comm, int * yac_instance_id )
  int yac_cget_default_instance_id()
  void yac_ccleanup_instance (int yac_instance_id)
  void yac_cdef_comp_instance ( int yac_instance_id,
                                const char * comp_name,
                                int * comp_id )
  void yac_cdef_comps_instance ( int yac_instance_id,
                                 const char ** comp_names,
                                 int num_comps,
                                 int * comp_ids )
  void yac_cpredef_comp_instance ( int yac_instance_id,
                                   const char * comp_name,
                                   int * comp_id )
  void yac_cget_comp_comm ( int comp_id, MPI_Comm* comp_comm )
  void yac_cdef_datetime_instance ( int yac_instance_id,
                                       const char * start_datetime,
                                       const char * end_datetime )
  void yac_cget_comps_comm_instance( int yac_instance_id,
                                     const char ** comp_names,
                                     int num_comps,
                                     MPI_Comm * comps_comm)
  void yac_cdef_calendar ( int calendar )
  int yac_cget_calendar ( )
  void yac_cenddef_instance ( int yac_instance_id )
  int yac_cget_nbr_comps_instance ( int yac_instance_id )
  int yac_cget_nbr_grids_instance ( int yac_instance_id )
  int yac_cget_comp_nbr_grids_instance ( int yac_instance_id, const char* comp_name )
  int yac_cget_nbr_fields_instance ( int yac_instance_id, const char * comp_name,
                                     const char* grid_name)
  void yac_cget_comp_names_instance ( int yac_instance_id, int nbr_comps,
                                      const char ** comp_names )
  void yac_cget_grid_names_instance ( int yac_instance_id, int nbr_grids,
                                      const char ** grid_names )
  void yac_cget_comp_grid_names_instance ( int yac_instance_id, const char* comp_name,
                                           int nbr_grids, const char ** grid_names )
  void yac_cget_field_names_instance ( int yac_instance_id, const char* comp_name,
                                       const char* grid_name,
                                       int nbr_fields, const char ** field_names )
  int yac_cget_field_id_instance ( int yac_instance_id, const char* comp_name,
                                   const char* grid_name,
                                   const char * field_name )
  const char* yac_cget_field_timestep_instance ( int yac_instance_id, const char* comp_name,
                                                 const char* grid_name,
                                                 const char * field_name )
  int yac_cget_field_role_instance ( int yac_instance_id, const char* comp_name,
                                     const char* grid_name, const char* field_name )
  void yac_cenable_field_frac_mask_instance ( int yac_instance_id,
                                              const char* comp_name,
                                              const char* grid_name,
                                              const char * field_name,
                                              double frac_mask_fallback_value)
  int yac_cget_field_collection_size_instance ( int yac_instance_id,
                                                const char* comp_name,
                                                const char* grid_name,
                                                const char * field_name )
  double yac_cget_field_frac_mask_fallback_value_instance ( int yac_instance_id,
                                                            const char* comp_name,
                                                            const char* grid_name,
                                                            const char * field_name )
  void yac_cdef_component_metadata_instance ( int yac_instance_id,
                                              const char* comp_name,
                                              const char* metadata)
  void yac_cdef_grid_metadata_instance ( int yac_instance_id,
                                         const char* grid_name,
                                         const char* metadata)
  void yac_cdef_field_metadata_instance ( int yac_instance_id,
                                          const char* comp_name,
                                          const char* grid_name,
                                          const char* field_name,
                                          const char* metadata)
  const char* yac_cget_component_metadata_instance(int yac_instance_id,
                                                   const char* comp_name)
  const char* yac_cget_grid_metadata_instance(int yac_instance_id,
                                              const char* grid_name)
  const char* yac_cget_field_metadata_instance(int yac_instance_id,
                                               const char* comp_name,
                                               const char* grid_name,
                                               const char* field_name)
  char * yac_cget_start_datetime_instance ( int yac_instance_id )
  char * yac_cget_end_datetime_instance ( int yac_instance_id )
  char * yac_cget_version ()
  void yac_cdef_grid_reg2d ( const char * grid_name,
                             int nbr_vertices[2],
                             int cyclic[2],
                             double *x_vertices,
                             double *y_vertices,
                             int *grid_id)
  void yac_cdef_points_reg2d ( const int grid_id,
                               int nbr_points[2],
                               const int location,
                               const double *x_points,
                               const double *y_points,
                               int *point_id )
  void yac_cdef_grid_curve2d ( const char * grid_name,
                               int nbr_vertices[2],
                               int cyclic[2],
                               double *x_vertices,
                               double *y_vertices,
                               int *grid_id)
  void yac_cdef_points_curve2d ( const int grid_id,
                                 int nbr_points[2],
                                 const int location,
                                 const double *x_points,
                                 const double *y_points,
                                 int *point_id )
  void yac_cdef_grid_cloud ( const char * grid_name,
                             int nbr_points,
                             double *x_points,
                             double *y_points,
                             int *grid_id)
  void yac_cdef_grid_unstruct ( const char * grid_name,
                                int nbr_vertices,
                                int nbr_cells,
                                int *num_vertices_per_cell,
                                double *x_vertices,
                                double *y_vertices,
                                int *cell_to_vertex,
                                int *grid_id)
  void yac_cdef_grid_unstruct_ll ( const char * grid_name,
                                   int nbr_vertices,
                                   int nbr_cells,
                                   int *num_vertices_per_cell,
                                   double *x_vertices,
                                   double *y_vertices,
                                   int *cell_to_vertex,
                                   int *grid_id)
  void yac_cdef_points_unstruct ( const int grid_id,
                                  const int nbr_points,
                                  const int location,
                                  const double *x_points,
                                  const double *y_points,
                                  int *point_id )
  void yac_cset_global_index ( const int * global_index,
                               int location,
                               int grid_id)
  void yac_cdef_field ( const char * field_name,
                        const int component_id,
                        const int * point_ids,
                        const int num_pointsets,
                        int collection_size,
                        const char* timestep,
                        int timeunit,
                        int * field_id )
  void yac_cdef_field_mask ( const char * field_name,
                             const int component_id,
                             const int * point_ids,
                             const int * mask_ids,
                             const int num_pointsets,
                             int collection_size,
                             const char* timestep,
                             int timeunit,
                             int * field_id )
  void yac_csync_def_instance ( int yac_instance_id )
  void yac_cget_ext_couple_config(int * ext_couple_config_id)
  void yac_cset_ext_couple_config_weight_file(int ext_couple_config_id,
                                              const char * weight_file)
  void yac_cset_ext_couple_config_mapping_side(int ext_couple_config_id,
                                               int mapping_side)
  void yac_cset_ext_couple_config_scale_factor(int ext_couple_config_id,
                                               double scale_factor)
  void yac_cset_ext_couple_config_scale_summand(int ext_couple_config_id,
                                                double scale_summand)
  void yac_cset_ext_couple_config_src_mask_names(
      int ext_couple_config_id, size_t num_src_mask_names,
      const char * const * src_mask_names)
  void yac_cset_ext_couple_config_tgt_mask_name(
      int ext_couple_config_id, const char * tgt_mask_name)
  void yac_cfree_ext_couple_config(int ext_couple_config_id)
  void yac_cdef_couple_custom_instance(int yac_instance_id,
                                       const char * src_comp_name,
                                       const char * src_grid_name,
                                       const char * src_field_name,
                                       const char * tgt_comp_name,
                                       const char * tgt_grid_name,
                                       const char * tgt_field_name,
                                       const char * coupling_timestep,
                                       int time_unit, int time_reduction,
                                       int interp_stack_config_id,
                                       int src_lag, int tgt_lag,
                                       int ext_couple_config_id)
  void yac_cget_ ( const int field_id,
                   const int collection_size,
                   double *recv_field,
                   int    *info,
                   int    *ierror ) nogil
  void yac_cget_async_ ( const int field_id,
                         const int collection_size,
                         double *recv_field,
                         int    *info,
                         int    *ierror )
  void yac_cput_ ( const int field_id,
                   const int collection_size,
                   double * send_field,
                   int *info,
                   int *ierror ) nogil
  void yac_cput_frac_ ( const int field_id,
                        const int collection_size,
                        double *send_field,
                        double *send_frac_mask,
                        int    *info,
                        int    *ierr ) nogil
  void yac_cexchange_ ( const int send_field_id,
                        const int recv_field_id,
                        const int collection_size,
                        double *send_field ,
                        double *recv_field,
                        int    *send_info,
                        int    *recv_info,
                        int    *ierror ) nogil
  void yac_cexchange_frac_ ( const int send_field_id,
                             const int recv_field_id,
                             const int collection_size,
                             double *send_field,
                             double *send_frac_mask,
                             double *recv_field,
                             int    *send_info,
                             int    *recv_info,
                             int    *ierror ) nogil
  void yac_ctest ( int field_id, int * flag )
  void yac_cwait ( int field_id ) nogil
  void yac_cupdate( int field_id )
  const char* yac_cget_field_name_from_field_id ( int field_id )
  const char* yac_cget_component_name_from_field_id ( int field_id )
  const char* yac_cget_grid_name_from_field_id ( int field_id )
  int yac_cget_role_from_field_id ( int field_id )
  const char* yac_cget_timestep_from_field_id ( int field_id )
  size_t yac_cget_grid_size ( int location, int grid_id )
  size_t yac_cget_points_size ( int points_id  )
  int yac_cget_collection_size_from_field_id ( const int field_id )
  const char* yac_cget_field_datetime(int field_id)
  void yac_cget_interp_stack_config(int * interp_stack_config_id)
  void yac_cadd_interp_stack_config_average(
  int interp_stack_config_id, int reduction_type, int partial_coverage)
  void yac_cadd_interp_stack_config_ncc(
  int interp_stack_config_id, int weight_type, int partial_coverage)
  void yac_cadd_interp_stack_config_nnn(int interp_stack_config_id, int type,
                                        unsigned int n,
                                        double max_search_distance,
                                        double scale)
  void yac_cadd_interp_stack_config_conservative(
      int interp_stack_config_id, int order, int enforced_conserv,
      int partial_coverage, int normalisation)
  void yac_cadd_interp_stack_config_spmap(
      int interp_stack_config_id, double spread_distance,
      double max_search_distance, int weight_type, int scale_type,
      double src_sphere_radius, double tgt_sphere_radius)
  void yac_cadd_interp_stack_config_hcsbb(int interp_stack_config_id)
  void yac_cadd_interp_stack_config_user_file(
      int interp_stack_config_id, char * filename)
  void yac_cadd_interp_stack_config_fixed(
      int interp_stack_config_id, double value)
  void yac_cadd_interp_stack_config_check(
      int interp_stack_config_id, char * constructor_key, char * do_search_key)
  void yac_cadd_interp_stack_config_creep(
      int interp_stack_config_id, int creep_distance)
  void yac_cget_interp_stack_config_from_string_yaml(char * interp_stack_config,
                                                     int * interp_stack_config_id)
  void yac_cget_interp_stack_config_from_string_json(char * interp_stack_config,
                                                     int * interp_stack_config_id)
  void yac_cfree_interp_stack_config(int interp_stack_config_id)
  void yac_cset_core_mask ( const int * is_core,
                            int location,
                            int grid_id)
  void yac_cset_mask ( const int * is_valid,
                       int points_id )
  void yac_cdef_mask_named ( const int grid_id,
                             const int nbr_points,
                             const int location,
                             const int * is_valid,
                             const char * name,
                             int *mask_id )
  void yac_cfinalize()
  ctypedef void (*yac_abort_func)(MPI_Comm comm, const char *msg,
                                  const char *source, int line) except *
  void yac_set_abort_handler(yac_abort_func custom_abort)
  yac_abort_func yac_get_abort_handler()
  void yac_cread_config_yaml_instance( int yac_instance_id,
                                       const char * yaml_file)
  void yac_cget_action( int field_id, int* action)
  void yac_cset_config_output_file_instance( int yac_instance_id,
                                             const char * filename,
                                             int fileformat,
                                             int sync_location,
                                             int include_definitions);
  void yac_ccompute_grid_cell_areas(int grid_id, double* cell_areas);

# helper functions for py interface (not in yac_interface.h)
cdef extern void yac_cget_comp_size_c2py(int comp_id, int* size)
cdef extern void yac_cget_comp_rank_c2py(int comp_id, int* rank)

_logger = logging.getLogger("yac")
_logger.addHandler(logging.NullHandler())

class Location(Enum):
    """
    Location for points

    Refers to @ref YAC_LOCATION_CELL, @ref YAC_LOCATION_CORNER and @ref YAC_LOCATION_EDGE
    """
    CELL = _LOCATION_CELL
    CORNER = _LOCATION_CORNER
    EDGE = _LOCATION_EDGE

class ExchangeType(Enum):
    """
    Exchange type of a field

    Refers to @ref YAC_EXCHANGE_TYPE_NONE, @ref YAC_EXCHANGE_TYPE_SOURCE and @ref YAC_EXCHANGE_TYPE_TARGET
    """
    NONE = _EXCHANGE_TYPE_NONE
    SOURCE = _EXCHANGE_TYPE_SOURCE
    TARGET = _EXCHANGE_TYPE_TARGET

class Action(Enum):
    """
    Refers to @ref YAC_ACTION_NONE, @ref YAC_ACTION_REDUCTION etc.
    """
    NONE = _ACTION_NONE
    REDUCTION = _ACTION_REDUCTION
    COUPLING = _ACTION_COUPLING
    GET_FOR_RESTART = _ACTION_GET_FOR_RESTART
    PUT_FOR_RESTART = _ACTION_PUT_FOR_RESTART
    OUT_OF_BOUND = _ACTION_OUT_OF_BOUND

class Reduction(Enum):
    """
    Reduction type for the definition of interpolations

    Refers to @ref YAC_REDUCTION_TIME_NONE, @ref YAC_REDUCTION_TIME_ACCUMULATE etc.
    """
    TIME_NONE = _REDUCTION_TIME_NONE
    TIME_ACCUMULATE = _REDUCTION_TIME_ACCUMULATE
    TIME_AVERAGE = _REDUCTION_TIME_AVERAGE
    TIME_MINIMUM = _REDUCTION_TIME_MINIMUM
    TIME_MAXIMUM = _REDUCTION_TIME_MAXIMUM

class Calendar(Enum):
    """
    Calendar type for use in def_calendar

    Refers to @ref YAC_CALENDAR_NOT_SET, @ref YAC_PROLEPTIC_GREGORIAN etc.
    """
    CALENDAR_NOT_SET = _CALENDAR_NOT_SET
    PROLEPTIC_GREGORIAN = _PROLEPTIC_GREGORIAN
    YEAR_OF_365_DAYS = _YEAR_OF_365_DAYS
    YEAR_OF_360_DAYS = _YEAR_OF_360_DAYS

def def_calendar(calendar : Calendar):
    """
    @see yac_cdef_calendar
    """
    yac_cdef_calendar(Calendar(calendar).value)

def get_calendar():
    """
    @see yac_cget_calendar
    """
    return Calendar(yac_cget_calendar())


class TimeUnit(Enum):
    """
    Refers to @ref YAC_TIME_UNIT_MILLISECOND, @ref YAC_TIME_UNIT_SECOND etc.
    """
    MILLISECOND = _TIME_UNIT_MILLISECOND
    SECOND = _TIME_UNIT_SECOND
    MINUTE = _TIME_UNIT_MINUTE
    HOUR = _TIME_UNIT_HOUR
    DAY = _TIME_UNIT_DAY
    MONTH = _TIME_UNIT_MONTH
    YEAR = _TIME_UNIT_YEAR
    ISO_FORMAT = _TIME_UNIT_ISO_FORMAT

class AverageReductionType(Enum):
    """
    Reduction type for average interpolation

    Refers to @ref YAC_AVG_ARITHMETIC, @ref YAC_AVG_DIST and @ref YAC_AVG_BARY
    """
    AVG_ARITHMETIC = _AVG_ARITHMETIC
    AVG_DIST = _AVG_DIST
    AVG_BARY = _AVG_BARY

class NCCReductionType(Enum):
    """
    Reduction type for ncc interpolation

    Refers to @ref YAC_NCC_AVG and @ref YAC_NCC_DIST
    """
    AVG = _NCC_AVG
    DIST = _NCC_DIST

class NNNReductionType(Enum):
    """
    Reduction type for nnn interpolation

    Refers to @ref YAC_NNN_AVG, @ref YAC_NNN_DIST etc.
    """
    AVG = _NNN_AVG
    DIST = _NNN_DIST
    GAUSS = _NNN_GAUSS
    RBF = _NNN_RBF

class ConservNormalizationType(Enum):
    """
    Normalization type for conservative interpolation

    Refers to @ref YAC_CONSERV_DESTAREA and @ref YAC_CONSERV_FRACAREA
    """
    DESTAREA = _CONSERV_DESTAREA
    FRACAREA = _CONSERV_FRACAREA

class SPMAPWeightType(Enum):
    """
    Refers to @ref YAC_SPMAP_AVG and @ref YAC_SPMAP_DIST
    """
    AVG = _SPMAP_AVG
    DIST = _SPMAP_DIST

class SPMAPScaleType(Enum):
    """
    Refers to @ref YAC_SPMAP_NONE, YAC_SPMAP_SRCAREA
    YAC_SPMAP_INVTGTAREA, and @ref YAC_SPMAP_FRACAREA
    """
    NONE = _SPMAP_NONE
    SRCAREA = _SPMAP_SRCAREA
    INVTGTAREA = _SPMAP_INVTGTAREA
    FRACAREA = _SPMAP_FRACAREA

class ConfigOutputFormat(Enum):
    YAML = _CONFIG_OUTPUT_FORMAT_YAML
    JSON = _CONFIG_OUTPUT_FORMAT_JSON

class ConfigOutputSyncLoc(Enum):
    DEF_COMP = _CONFIG_OUTPUT_SYNC_LOC_DEF_COMP # after component definition
    SYNC_DEF = _CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF # after synchronization of definition
    ENDDEF   = _CONFIG_OUTPUT_SYNC_LOC_ENDDEF # after end of definitions

class YAC:
    """
    Initializies a YAC instance and provides further functionality

    The destructor finalizes the YAC instance by calling yac_cfinalize_instance
    """
    def __init__(self, comm = None, default_instance = False):
        """
        @see yac_cinit_instance
        """
        cdef int instance_id
        cdef MPI_Comm c_comm
        if comm is None:
            if default_instance:
                _logger.debug("init")
                yac_cinit()
                instance_id = yac_cget_default_instance_id()
            else:
                _logger.debug("init_instance")
                yac_cinit_instance(&instance_id)
        else:
            from mpi4py import MPI
            if type(comm) is MPI.Intracomm:
                comm = MPI.Comm.py2f(comm)
            if default_instance:
                _logger.debug("init_comm")
                yac_cinit_comm(MPI_Comm_f2c(comm))
                instance_id = yac_cget_default_instance_id()
            else:
                _logger.debug("init_comm_instance")
                yac_cinit_comm_instance(MPI_Comm_f2c(comm), &instance_id)
        _logger.debug(f"instance_id={instance_id}")
        self.instance_id = instance_id
        self.owned_instance = True

    @classmethod
    @property
    def default_instance(cls):
        return cls.from_id(yac_cget_default_instance_id())

    @classmethod
    def from_id(cls, id):
        yac = cls.__new__(cls)
        yac.instance_id = id
        yac.owned_instance = False
        return yac

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        """
        @see yac_ccleanup_instance
        """
        if self.owned_instance:
            _logger.debug(f"cleanup instance_id={self.instance_id}")
            yac_ccleanup_instance(self.instance_id)
            self.owned_instance = False

    def def_comp(self, comp_name : str):
        """
        @see yac_cdef_comp_instance
        """
        cdef int comp_id
        _logger.debug(f"def_comp: comp_name={comp_name}")
        yac_cdef_comp_instance(self.instance_id, comp_name.encode(), &comp_id)
        _logger.debug(f"comp_id={comp_id}")
        return Component(comp_id)

    def def_comps(self, comp_names = []):
        """
        @see yac_cdef_comps_instance
        """
        cdef int comp_len = len(comp_names)
        cdef const char **c_comp_names = <const char **>malloc(comp_len * sizeof(const char *))
        cdef int *c_comp_ids = <int*>malloc(comp_len * sizeof(int))
        _logger.debug(f"def_comps: comp_names={comp_names}")
        byte_comp_names = [c.encode() for c in comp_names]
        for i in range(comp_len):
            c_comp_names[i] = byte_comp_names[i]
        yac_cdef_comps_instance(self.instance_id, c_comp_names, comp_len, c_comp_ids)
        comp_list = [Component(c_comp_ids[i]) for i in range(comp_len) ]
        free(c_comp_names)
        free(c_comp_ids)
        _logger.debug("comp_list={comp_list}")
        return comp_list

    def predef_comp(self, comp_name :str):
        """
        @see yac_cpredef_comp_instance
        """
        cdef int comp_id
        _logger.debug(f"predef_comp: comp_name={comp_name}")
        yac_cpredef_comp_instance(self.instance_id, comp_name.encode(), &comp_id)
        _logger.debug("comp_id={comp_id}")
        return Component(comp_id)

    def def_datetime(self, start_datetime, end_datetime):
        """
        @see yac_cdef_datetime_instance

        The parameters can be given either as a string in iso8601
        format or as datetime objects
        """
        try:
            import datetime
            if(type(start_datetime) is datetime.datetime):
                start_datetime = start_datetime.isoformat()
        except:
            pass
        if start_datetime is not None:
            _logger.debug(f"def_datetime start: {start_datetime}")
            yac_cdef_datetime_instance ( self.instance_id,
                                         start_datetime.encode(),
                                         NULL)

        try:
            if(type(end_datetime) is datetime.datetime):
                end_datetime = end_datetime.isoformat()
        except:
            pass
        if end_datetime is not None:
            _logger.debug(f"def_datetime end: {end_datetime}")
            yac_cdef_datetime_instance ( self.instance_id,
                                         NULL,
                                         end_datetime.encode())

    @property
    def start_datetime(self):
        """
        @see yac_cget_start_datetime_instance (`datetime.datetime`, read-only).
        """
        start = yac_cget_start_datetime_instance(self.instance_id)
        return bytes.decode(start)

    @property
    def end_datetime(self):
        """
        @see yac_cget_end_datetime_instance (`datetime.datetime`, read-only).
        """
        end = yac_cget_end_datetime_instance(self.instance_id)
        return bytes.decode(end)

    def sync_def(self):
        """
        @see yac_csync_def_instance
        """
        _logger.debug("sync_def")
        yac_csync_def_instance(self.instance_id)

    def def_couple(self,
                   src_comp : str, src_grid : str, src_field,
                   tgt_comp : str, tgt_grid : str, tgt_field,
                   coupling_timestep : str, timeunit : TimeUnit,
                   time_reduction : Reduction,
                   interp_stack, src_lag = 0, tgt_lag = 0,
                   weight_file = None, mapping_on_source = 1,
                   scale_factor = 1.0, scale_summand = 0.0,
                   src_masks_names = None, tgt_mask_name = None):
        """
        @see yac_cdef_couple_instance
        """
        cdef char * weight_file_ptr
        if weight_file is None:
            weight_file_ptr = NULL
        else:
            weight_file_bytes = weight_file.encode()
            weight_file_ptr = weight_file_bytes
        cdef const char ** src_mask_names_ptr = NULL
        cdef const char * tgt_mask_name_ptr = NULL
        if tgt_mask_name is not None:
            tgt_mask_enc = tgt_mask_name.encode()
            tgt_mask_name = tgt_mask_enc
        cdef int couple_config_id
        yac_cget_ext_couple_config(&couple_config_id)
        yac_cset_ext_couple_config_weight_file(couple_config_id,
                                               weight_file_ptr)
        yac_cset_ext_couple_config_mapping_side(couple_config_id,
                                                mapping_on_source)
        yac_cset_ext_couple_config_scale_factor(couple_config_id,
                                                scale_factor)
        yac_cset_ext_couple_config_scale_summand(couple_config_id,
                                                 scale_summand)
        if src_masks_names is not None:
            if type(src_masks_names) is str:
                src_masks = [src_masks_names]
            src_masks_enc = [s.encode() for s in src_masks_names]
            src_mask_names_ptr = <const char **>malloc(len(src_masks_enc) * sizeof(char*))
            for i in range(len(src_masks_enc)):
                src_mask_names_ptr[i] = src_masks_enc[i]
            yac_cset_ext_couple_config_src_mask_names(couple_config_id,
                                                      len(src_masks_enc),
                                                      src_mask_names_ptr)
            free(src_mask_names_ptr)
        if tgt_mask_name is not None:
            yac_cset_ext_couple_config_tgt_mask_name(couple_config_id, tgt_mask_name.encode())
        _logger.debug(f"def_couple {(src_comp.encode(), src_grid.encode(), src_field.encode())}, {(tgt_comp.encode(), tgt_grid.encode(), tgt_field.encode())}")
        yac_cdef_couple_custom_instance(self.instance_id,
                                        src_comp.encode(), src_grid.encode(), src_field.encode(),
                                        tgt_comp.encode(), tgt_grid.encode(), tgt_field.encode(),
                                        coupling_timestep.encode(), TimeUnit(timeunit).value,
                                        Reduction(time_reduction).value,
                                        interp_stack.interp_stack_id, src_lag, tgt_lag,
                                        couple_config_id)
        yac_cfree_ext_couple_config(couple_config_id)

    def enddef(self):
        """
        @see yac_cenddef_instance
        """
        _logger.debug("enddef")
        yac_cenddef_instance(self.instance_id)

    @property
    def component_names(self):
        """
        @see yac_cget_comp_names
        """
        cdef int nbr_components = yac_cget_nbr_comps_instance(self.instance_id)
        cdef const char **ret = <const char **>malloc(nbr_components * sizeof(const char *))
        yac_cget_comp_names_instance(self.instance_id, nbr_components, ret)
        comp_list = [bytes(ret[i]).decode('UTF-8') for i in range(nbr_components) ]
        free(ret)
        return comp_list

    @property
    def grid_names(self):
        """
        @see yac_cget_grid_names
        """
        cdef int nbr_grids = yac_cget_nbr_grids_instance(self.instance_id)
        cdef const char **ret = <const char **>malloc(nbr_grids * sizeof(const char *))
        yac_cget_grid_names_instance(self.instance_id, nbr_grids, ret)
        grid_list = [bytes(ret[i]).decode('UTF-8') for i in range(nbr_grids) ]
        free(ret)
        return grid_list

    def get_comp_grid_names(self, comp_name):
        """
        @see yac_cget_comp_grid_names
        """
        cdef int nbr_grids = yac_cget_comp_nbr_grids_instance(self.instance_id, comp_name.encode())
        cdef const char **ret = <const char **>malloc(nbr_grids * sizeof(const char *))
        yac_cget_comp_grid_names_instance(self.instance_id, comp_name.encode(), nbr_grids, ret)
        grid_list = [bytes(ret[i]).decode('UTF-8') for i in range(nbr_grids) ]
        free(ret)
        return grid_list

    def get_field_names(self, comp_name : str, grid_name : str):
        """
        @see yac_cget_field_names
        """
        cdef int nbr_fields = yac_cget_nbr_fields_instance(self.instance_id,
                                                                 comp_name.encode(),
                                                                 grid_name.encode())
        cdef const char **ret = <const char **>malloc(nbr_fields * sizeof(const char *))
        yac_cget_field_names_instance(self.instance_id, comp_name.encode(),
                                            grid_name.encode(), nbr_fields, ret)
        field_list = [bytes(ret[i]).decode('UTF-8') for i in range(nbr_fields) ]
        free(ret)
        return field_list

    def get_field_id(self, comp_name : str, grid_name : str, field_name : str):
        """
        @see yac_cget_field_id
        """
        return yac_cget_field_id_instance (self.instance_id,
                                           comp_name.encode(),
                                           grid_name.encode(),
                                           field_name.encode())

    def get_field_timestep(self, comp_name : str, grid_name : str, field_name : str):
        """
        @see yac_cget_field_timestep
        """
        return yac_cget_field_timestep_instance(self.instance_id,
                                                comp_name.encode(),
                                                grid_name.encode(),
                                                field_name.encode()).decode('UTF-8')

    def get_field_role(self, comp_name : str, grid_name : str, field_name : str):
        """
        @see yac_cget_field_role
        """
        return ExchangeType(yac_cget_field_role_instance (self.instance_id,
                                                          comp_name.encode(),
                                                          grid_name.encode(),
                                                          field_name.encode()))

    def get_field_collection_size(self, comp_name : str, grid_name : str, field_name : str):
        """
        @see yac_cget_field_collection_size
        """
        return yac_cget_field_collection_size_instance(self.instance_id,
                                                       comp_name.encode(),
                                                       grid_name.encode(),
                                                       field_name.encode())

    def get_field_frac_mask_fallback_value(self, comp_name : str, grid_name : str, field_name : str):
        """
        @see yac_cget_field_frac_mask_fallback_value
        """
        return yac_cget_field_frac_mask_fallback_value_instance(self.instance_id,
                                                                comp_name.encode(),
                                                                grid_name.encode(),
                                                                field_name.encode())

    def enable_field_frac_mask(self, comp_name : str, grid_name : str, field_name : str,
                               frac_mask_fallback_value : _np.float64):
        """
        @see yac_cenable_field_frac_mask
        """
        _logger.debug(f"enable_field_frac_mask {(comp_name, grid_name, field_name)})")
        yac_cenable_field_frac_mask_instance(self.instance_id,
                                             comp_name.encode(),
                                             grid_name.encode(),
                                             field_name.encode(),
                                             frac_mask_fallback_value)

    def def_component_metadata(self, comp_name : str, metadata : bytes):
        """
        @see yac_cdef_component_metadata
        """
        _logger.debug(f"def_component_metadata comp_name={comp_name}")
        yac_cdef_component_metadata_instance(self.instance_id ,
                                             comp_name.encode(), metadata)

    def def_grid_metadata(self, grid_name : str, metadata : bytes):
        """
        @see yac_cdef_grid_metadata
        """
        _logger.debug(f"def_grid_metadata grid_name={grid_name}")
        yac_cdef_grid_metadata_instance(self.instance_id,
                                        grid_name.encode(), metadata)

    def def_field_metadata(self, comp_name : str, grid_name : str,
                           field_name : str,metadata : bytes):
        """
        @see yac_cdef_field_metadata
        """
        _logger.debug("def_field_metadata (comp_name, grid_name, field_name)="
                      f"{(comp_name, grid_name, field_name)}")
        yac_cdef_field_metadata_instance(self.instance_id, comp_name.encode(),
                                         grid_name.encode(), field_name.encode(),
                                         metadata)

    def get_component_metadata(self, comp_name : str):
        """
        @see yac_cget_component_metadata
        """
        cdef const char* metadata = yac_cget_component_metadata_instance(self.instance_id,
                                                                         comp_name.encode())
        return bytes(metadata).decode('UTF-8') if metadata != NULL else None

    def get_grid_metadata(self, grid_name : str):
        """
        @see yac_cget_grid_metadata
        """
        cdef const char* metadata = yac_cget_grid_metadata_instance(self.instance_id,
                                                                    grid_name.encode())
        return bytes(metadata).decode('UTF-8') if metadata != NULL else None

    def get_field_metadata(self, comp_name : str, grid_name : str, field_name :str):
        """
        @see yac_cget_field_metadata
        """
        cdef const char* metadata = yac_cget_field_metadata_instance(self.instance_id,
                                                                     comp_name.encode(),
                                                                     grid_name.encode(),
                                                                     field_name.encode())
        return bytes(metadata).decode('UTF-8') if metadata != NULL else None

    def get_comps_comm(self, comp_names):
        """
        @see yac_cget_comps_comm
        """
        from mpi4py import MPI
        cdef MPI_Comm comm
        cptr = [c.encode() for c in comp_names]
        cdef const char ** comp_names_c_ptr = <const char **>malloc(len(comp_names) * sizeof(const char *))
        for i in range(len(comp_names)):
            comp_names_c_ptr[i] = cptr[i]
        yac_cget_comps_comm_instance(self.instance_id, comp_names_c_ptr, len(comp_names), &comm)
        free(comp_names_c_ptr)
        # convert to mpi4py communicator
        return MPI.Comm.f2py(MPI_Comm_c2f(comm))

    def read_config_yaml(self, yaml_file : str):
        """
        @see yac_cread_config_yaml_instance
        """
        _logger.debug(f"read_config_yaml yaml_file={yaml_file}")
        yac_cread_config_yaml_instance(self.instance_id, yaml_file.encode())

    def set_config_output_file(self, filename : str,
                               fileformat : ConfigOutputFormat,
                               sync_location : ConfigOutputSyncLoc,
                               include_definitions : bool = False):
        """
        @see yac_cset_config_output_file_instance
        """
        yac_cset_config_output_file_instance(self.instance_id,
                                             filename.encode(),
                                             fileformat.value,
                                             sync_location.value,
                                             1 if include_definitions else 0)

class Component:
    """
    Stores the component_id and provides further functionality
    """
    def __init__(self, comp_id):
        self.comp_id = comp_id

    @property
    def comp_comm(self):
        """
        @see yac_cget_comp_comm (`MPI.Comm`, read-only)
        """
        from mpi4py import MPI
        cdef MPI_Comm comm
        yac_cget_comp_comm(self.comp_id, &comm)
        # convert to mpi4py communicator
        return MPI.Comm.f2py(MPI_Comm_c2f(comm))

    @property
    def size(self):
        """
        number of processes in this component
        """
        cdef int size
        yac_cget_comp_size_c2py(self.comp_id, &size)
        return size

    @property
    def rank(self):
        """
        process index in the component
        """
        cdef int rank
        yac_cget_comp_rank_c2py(self.comp_id, &rank)
        return rank

class Mask:
    """
    Stores the mask_id
    """
    def __init__(self, mask_id):
        self.mask_id = mask_id

class Grid:
    """
    Stores the grid_id and provides further functionality

    Base class for Reg2dGrid and UnstructuredGrid
    """
    def __init__(self, grid_id):
        self.grid_id = grid_id

    def set_global_index(self, global_index, location : Location):
        """
        @see yac_cset_global_index
        """
        cdef const int[::1] global_index_view = _np.ascontiguousarray(global_index, dtype=_np.intc)
        assert len(global_index_view) == yac_cget_grid_size ( location.value, self.grid_id ), "Wrong number of indices provided"
        _logger.debug(f"set_global_index grid_id={self.grid_id}")
        yac_cset_global_index(&global_index_view[0], location.value, self.grid_id)

    @property
    def nbr_cells(self):
        """
        @see yac_cget_grid_size (`int`, read-only)
        """
        return yac_cget_grid_size ( Location.CELL.value, self.grid_id )

    @property
    def nbr_corners(self):
        """
        @see yac_cget_grid_size (`int`, read-only)
        """
        return yac_cget_grid_size ( Location.CORNER.value, self.grid_id )

    @property
    def nbr_edges(self):
        """
        @see yac_cget_grid_size (`int`, read-only)
        """
        return yac_cget_grid_size ( Location.EDGE.value, self.grid_id )

    def set_core_mask(self, is_core, location : Location):
        """
        @see yac_cset_core_mask
        """
        cdef size_t len_is_core = len(is_core)
        assert(len_is_core == yac_cget_grid_size( location.value, self.grid_id ) )
        cdef const int[::1] np_mask = _np.ascontiguousarray(is_core, dtype=_np.intc)
        _logger.debug(f"set_core_mask grid_id={self.grid_id}")
        yac_cset_core_mask ( &np_mask[0], location.value, self.grid_id)

    def def_mask(self, location : Location,
                 is_valid, name = None):
        cdef int len_is_valid = len(is_valid)
        cdef const int[::1] np_mask = _np.ascontiguousarray(is_valid, dtype=_np.int32)
        cdef int mask_id
        cdef char* c_name = NULL
        if name is not None:
            name_enc = name.encode()
            c_name = name_enc
        _logger.debug(f"def_mask_named grid_id={self.grid_id}, name={name}")
        yac_cdef_mask_named ( self.grid_id,
                              len_is_valid,
                              location.value,
                              &np_mask[0],
                              c_name,
                              &mask_id )
        return Mask(mask_id)

    def compute_grid_cell_areas(self, cell_areas = None):
        """
        @see yac_ccompute_grid_cell_areas
        """
        if cell_areas is None:
            cell_areas = _np.empty(self.nbr_cells)
        cdef double[::1] np_cell_areas = _np.ascontiguousarray(cell_areas, dtype=_np.double)
        yac_ccompute_grid_cell_areas(self.grid_id, &np_cell_areas[0])
        return cell_areas

class Points:
    """
    Stores the points_id and provides further functionality
    """
    def __init__(self, points_id):
        self.points_id = points_id

    @property
    def size(self):
        """
        @see yac_cget_points_size (`int`, read-only)
        """
        return yac_cget_points_size ( self.points_id )

    def set_mask(self, is_valid):
        """
        @see yac_cset_mask
        """
        cdef size_t len_is_valid = len(is_valid)
        assert len_is_valid==self.size
        cdef const int[::1] np_mask = _np.ascontiguousarray(is_valid, dtype=_np.intc)
        _logger.debug(f"set_mask points_id={self.points_id}")
        yac_cset_mask ( &np_mask[0],
                        self.points_id )

class Reg2dGrid(Grid):
    """
    A stuctured 2d Grid
    """
    def __init__(self, grid_name : str, x_vertices, y_vertices,
                 cyclic = [False, False]):
        """
        @see yac_cdef_grid_reg2d
        """
        cdef int grid_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        cdef int[2] cyclic_view = cyclic
        _logger.debug(f"def_grid_reg2d grid_name={grid_name}")
        yac_cdef_grid_reg2d(grid_name.encode(), [len(x),len(y)],
                            cyclic_view, <double*>&x[0], <double*>&y[0], &grid_id)
        _logger.debug(f"grid_id={grid_id}")
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_reg2d
        """
        cdef int points_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        _logger.debug(f"def_points_reg2d grid_id={self.grid_id}")
        yac_cdef_points_reg2d(self.grid_id, [len(x), len(y)],
                              location.value, &x[0], &y[0], &points_id)
        _logger.debug(f"points_id={points_id}")
        return Points(points_id)

    def def_points_unstruct(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_unstruct
        """
        assert len(x_vertices) == len(y_vertices)
        cdef int points_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        _logger.debug(f"def_points_unstruct grid_id={self.grid_id}")
        yac_cdef_points_unstruct(self.grid_id, len(x),
                                 location.value, &x[0], &y[0], &points_id)
        _logger.debug(f"points_id={points_id}")
        return Points(points_id)

class Curve2dGrid(Grid):
    """
    A curvilinear stuctured 2d Grid
    """
    def __init__(self, grid_name : str, x_vertices, y_vertices,
                 cyclic = [False, False]):
        """
        @see yac_cdef_grid_curve2d
        """
        cdef int grid_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices.flatten(), dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices.flatten(), dtype=_np.double)
        cdef int[2] cyclic_view = cyclic
        _logger.debug(f"def_grid_curve2d grid_name={grid_name}")
        yac_cdef_grid_curve2d(grid_name.encode(),
                              [_np.shape(x_vertices)[1], _np.shape(y_vertices)[0]],
                              cyclic_view, <double*>&x[0], <double*>&y[0], &grid_id)
        _logger.debug(f"grid_id={grid_id}")
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_curve2d
        """
        assert x_vertices.shape == y_vertices.shape
        cdef int points_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices.flatten(), dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices.flatten(), dtype=_np.double)
        _logger.debug(f"def_points_curve2d grid_id={self.grid_id}")
        yac_cdef_points_curve2d(self.grid_id,
                                [_np.shape(x_vertices)[1], _np.shape(x_vertices)[0]],
                                location.value, &x[0], &y[0], &points_id)
        _logger.debug(f"points_id={points_id}")
        return Points(points_id)

    def def_points_unstruct(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_unstruct
        """
        assert len(x_vertices) == len(y_vertices)
        cdef int points_id
        cdef const double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef const double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        _logger.debug(f"def_points_unstruct grid_id={self.grid_id}")
        yac_cdef_points_unstruct(self.grid_id, len(x),
                                 location.value, &x[0], &y[0], &points_id)
        _logger.debug(f"points_id={points_id}")
        return Points(points_id)

class UnstructuredGrid(Grid):
    """
    An unstuctured 2d Grid
    """
    def __init__(self, grid_name : str, num_vertices_per_cell,
                 x_vertices, y_vertices, cell_to_vertex, use_ll_edges=False):
        """
        @see yac_cdef_grid_unstruct and @see yac_cdef_grid_unstruct_ll
        """
        cdef int grid_id
        cdef const int[::1] num_vertices_per_cell_view = _np.ascontiguousarray(num_vertices_per_cell, dtype=_np.intc)
        cdef const double[::1] x_vertices_view = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef const double[::1] y_vertices_view = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        cdef const int[::1] cell_to_vertex_view = _np.ascontiguousarray(cell_to_vertex, dtype=_np.intc)
        if len(num_vertices_per_cell_view) == 0:
            _logger.debug(f"def_grid_cloud grid_name={grid_name}")
            yac_cdef_grid_cloud (grid_name.encode(),
                                 len(x_vertices_view),
                                 <double*>&x_vertices_view[0],
                                 <double*>&y_vertices_view[0],
                                 &grid_id)
        elif not use_ll_edges:
            _logger.debug(f"def_grid_unstruct grid_name={grid_name}")
            yac_cdef_grid_unstruct(grid_name.encode(), len(x_vertices_view),
                                   len(num_vertices_per_cell_view),
                                   <int*>&num_vertices_per_cell_view[0],
                                   <double*>&x_vertices_view[0],
                                   <double*>&y_vertices_view[0],
                                   <int*>&cell_to_vertex_view[0], &grid_id)
        else:
            _logger.debug(f"def_grid_unstruct_ll grid_name={grid_name}")
            yac_cdef_grid_unstruct_ll(grid_name.encode(), len(x_vertices_view),
                                      len(num_vertices_per_cell_view),
                                      <int*>&num_vertices_per_cell_view[0],
                                      <double*>&x_vertices_view[0],
                                      <double*>&y_vertices_view[0],
                                      <int*>&cell_to_vertex_view[0], &grid_id)
        _logger.debug(f"grid_id={grid_id}")
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_points, y_points):
        """
        @see yac_cdef_points_unstruct
        """
        cdef int points_id
        cdef const double[::1] x_points_view = _np.ascontiguousarray(x_points, dtype=_np.double)
        cdef const double[::1] y_points_view = _np.ascontiguousarray(y_points, dtype=_np.double)
        _logger.debug(f"yac_cdef_points_unstruct grid_id={self.grid_id}")
        yac_cdef_points_unstruct(self.grid_id, len(x_points_view), location.value,
                                 &x_points_view[0], &y_points_view[0], &points_id)
        _logger.debug(f"points_id={points_id}")
        return Points(points_id)

class CloudGrid(UnstructuredGrid):
    """
    @see yac_cdef_grid_cloud
    """
    def __init__(self, grid_name : str, x_vertices, y_vertices):
        super().__init__(grid_name, [], x_vertices, y_vertices, [])

    def def_points(self, x_points, y_points):
        """
        @see yac_cdef_points_unstruct
        Uses YAC_LOCATION_CORNER as location
        """
        return super().def_points(Location.CORNER, x_points, y_points)

class Field:
    """
    Store the field_id
    """
    def __init__(self, field_id, size=None):
        self.field_id = field_id
        self._size = size

    @classmethod
    def create(cls, field_name : str, comp : Component, points, collection_size,
               timestep : str, timeunit : TimeUnit, masks = None):
        """
        @see yac_cdef_field
        """
        from collections.abc import Iterable
        cdef int field_id
        if not isinstance(points, Iterable):
            points = [points]
        cdef const int[:] point_ids_array = _np.array([p.points_id for p in points], dtype=_np.intc)
        size = sum(p.size for p in points)
        cdef const int[:] mask_ids_array
        if masks is None:
            _logger.debug(f"def_field field_name={field_name}, comp_id={comp.comp_id}, point_ids={[p.points_id for p in points]}")
            yac_cdef_field(field_name.encode(), comp.comp_id,
                           &point_ids_array[0], len(point_ids_array),
                           collection_size, timestep.encode(), TimeUnit(timeunit).value, &field_id)
        else:
            if not isinstance(masks, Iterable):
                masks = [masks]
            mask_ids_array = _np.array([m.mask_id for m in masks], dtype=_np.intc)
            _logger.debug(f"def_field_mask field_name={field_name}, comp_id={comp.comp_id}, point_ids={[p.points_id for p in points]}")
            yac_cdef_field_mask(field_name.encode(), comp.comp_id,
                                &point_ids_array[0], &mask_ids_array[0],
                                len(point_ids_array),
                                collection_size, timestep.encode(), TimeUnit(timeunit).value, &field_id)
        _logger.debug(f"field_id={field_id}")
        return Field(field_id, size)

    def test(self):
        """
        @see yac_ctest

        """
        cdef int flag
        _logger.debug(f"test field_id={self.field_id}")
        yac_ctest ( self.field_id, &flag )
        _logger.debug(f"flag={flag}")
        return flag

    def wait(self):
        """
        @see yac_cwait

        """
        cdef int field_id = self.field_id
        _logger.debug(f"wait field_id={field_id}")
        with cython.nogil:
            yac_cwait ( field_id )

    @coroutine
    def wait_coro(self):
        """
        Coroutine. Blocks until the communication is completed.
        """
        while self.test() == 0:
            yield

    @cython.boundscheck(False)
    def get(self, buf=None, asyn=False):
        """
        @see yac_cget_

        @param[out] buf    receive buffer, if `None` a numpy array of correct size is allocated
        @param[in] asyn    if True the call returns immidiatly and must be completed
                           with `test` or `wait`
        """
        cdef int info
        cdef int ierror
        buf_in = buf
        if buf is None:
            buf = _np.empty((self.collection_size, self.size), dtype=_np.double)
        buf = _np.ascontiguousarray(buf.reshape(self.collection_size, self.size), dtype=_np.double)
        if buf_in is not None and buf.base is None:
            _logger.warning("get: non-contiguous buffer passed to get. Reallocated memory.")
        cdef double[:,::1] buf_view = buf
        cdef int field_id = self.field_id
        cdef int collection_size = self.collection_size
        if asyn:
            _logger.debug(f"get_async_ field_id={self.field_id}")
            yac_cget_async_(field_id, collection_size, &buf_view[0,0], &info, &ierror)
        else:
            _logger.debug(f"get_ field_id={self.field_id}")
            with cython.nogil:
                yac_cget_(field_id, collection_size, &buf_view[0,0], &info, &ierror)
        if ierror != 0:
            raise RuntimeError("yac_cget returned error number " + str(ierror))
        _logger.debug(f"info={info}")
        return buf, Action(info)

    async def get_coro(self, buf=None):
        """
        Coroutine. Executes a get operation.

        @see yac_cget_async_

        @param[out] buf    receive buffer, if `None` a numpy array of correct size is allocated
        """
        buf, info = self.get(buf, asyn=True)
        await self.wait_coro()
        return buf, info

    @cython.boundscheck(False)
    def put(self, buf, frac_mask = None):
        """
        @see yac_cput_
        """
        cdef int info
        cdef int ierror
        cdef const double[:,:,::1] buf_view = _np.ascontiguousarray(buf.reshape(-1,self.collection_size,self.size), dtype=_np.double)
        cdef const double[:,:,::1] frac_mask_view
        cdef int field_id = self.field_id
        cdef int collection_size = self.collection_size
        if frac_mask is not None:
            frac_mask_view = _np.ascontiguousarray(
                frac_mask.reshape(-1,self.collection_size,self.size), dtype=_np.double)
            _logger.debug(f"put_frac_ field_id={self.field_id}")
            with cython.nogil:
                yac_cput_frac_(field_id, collection_size, <double*>&buf_view[0,0,0],
                               <double*>&frac_mask_view[0,0,0], &info, &ierror)
        else:
            _logger.debug(f"put_ field_id={self.field_id}")
            with cython.nogil:
                yac_cput_(field_id, collection_size, <double*>&buf_view[0,0,0], &info, &ierror)
        if ierror != 0:
            raise RuntimeError("yac_cput returned error number " + str(ierror))
        _logger.debug(f"info={info}")
        return Action(info)

    async def put_coro(self, buf, frac_mask = None):
        """
        Coroutine. Executes a put operation.

        @see yac_cput_
        """
        await self.wait_coro()
        return self.put(buf, frac_mask)

    def update(self):
        """
        @see yac_cupdate
        """
        _logger.debug(f"update field_id={self.field_id}")
        yac_cupdate(self.field_id)

    @classmethod
    @cython.boundscheck(False)
    def exchange(cls, send_field, recv_field, send_buf, recv_buf=None, send_frac_mask=None):
        """
        @see yac_cexchange_frac_
        """
        assert send_field.collection_size == recv_field.collection_size, "exchange can only be used with fields of the same collection_size"
        cdef int collection_size = send_field.collection_size
        cdef int send_info
        cdef int recv_info
        cdef int ierror
        cdef const double[:,:,::1] send_buf_view = _np.ascontiguousarray(send_buf.reshape(-1,collection_size,send_field.size), dtype=_np.double)
        cdef const double[:,:,::1] frac_mask_view
        if recv_buf is None:
            recv_buf = _np.empty((collection_size, recv_field.size), dtype=_np.double)
        recv_buf = _np.ascontiguousarray(recv_buf.reshape(collection_size, recv_field.size), dtype=_np.double)
        cdef double[:,::1] recv_buf_view = recv_buf
        cdef int send_field_id = send_field.field_id
        cdef int recv_field_id = recv_field.field_id
        if send_frac_mask is None:
            _logger.debug(f"exchange_ with {send_field.field_id=} and {recv_field.field_id=}")
            with cython.nogil:
                yac_cexchange_(send_field_id, recv_field_id, collection_size,
                               <double*>&send_buf_view[0,0,0],
                               <double*>&recv_buf_view[0,0],
                               &send_info, &recv_info, &ierror)
        else:
            _logger.debug(f"exchange_frac_ with {send_field.field_id=} and {recv_field.field_id=}")
            frac_mask_view = _np.ascontiguousarray(
                send_frac_mask.reshape(-1,collection_size,send_field.size), dtype=_np.double)
            with cython.nogil:
                yac_cexchange_frac_(send_field_id, recv_field_id, collection_size,
                                    <double*>&send_buf_view[0,0,0],
                                    <double*>&frac_mask_view[0,0,0],
                                    <double*>&recv_buf_view[0,0],
                                    &send_info, &recv_info, &ierror)
        if ierror != 0:
            raise RuntimeError("yac_cexchange_frac_ returned error number " + str(ierror))

        return recv_buf, send_info, recv_info

    @property
    def name(self):
        """
        @see yac_cget_field_name_from_field_id
        """
        return bytes.decode(yac_cget_field_name_from_field_id ( self.field_id ))

    @property
    def grid_name(self):
        """
        @see yac_cget_grid_name_from_field_id
        """
        return bytes.decode(yac_cget_grid_name_from_field_id ( self.field_id ))

    @property
    def component_name(self):
        """
        @see yac_cget_component_name_from_field_id
        """
        return bytes.decode(yac_cget_component_name_from_field_id ( self.field_id ))

    @property
    def role(self):
        """
        @see yac_cget_role_from_field_id
        """
        return ExchangeType(yac_cget_role_from_field_id ( self.field_id ))

    @property
    def timestep(self):
        """
        @see yac_cget_timestep_from_field_id
        """
        return bytes.decode(yac_cget_timestep_from_field_id ( self.field_id ))

    @property
    def collection_size(self):
        """
        @see yac_cget_collection_size_from_field_id
        """
        return yac_cget_collection_size_from_field_id(self.field_id)

    @property
    def size(self):
        """
        The size of the corresponding points object
        """
        return self._size

    @property
    def datetime(self):
        """
        @see yac_cget_field_datetime
        """
        return bytes.decode(yac_cget_field_datetime(self.field_id))

    @property
    def action(self):
        """
        @see yac_cget_action
        """
        cdef int action
        yac_cget_action(self.field_id, &action)
        return Action(action)

class InterpolationStack:
    def __init__(self):
        """
        @see yac_cget_interp_stack_config
        """
        cdef int interp_stack_config_id
        _logger.debug(f"get_interp_stack_config")
        yac_cget_interp_stack_config(&interp_stack_config_id)
        _logger.debug(f"interp_stack_config_id={interp_stack_config_id}")
        self.interp_stack_id = interp_stack_config_id

    @classmethod
    def from_string_yaml(cls, interp_stack_string : str):
        interp_stack = cls.__new__(cls)
        cdef int interp_stack_config_id
        _logger.debug(f"get_interp_stack_config_from_string_yaml: interp_stack_string={interp_stack_string}")
        yac_cget_interp_stack_config_from_string_yaml(interp_stack_string.encode(), &interp_stack_config_id)
        _logger.debug(f"interp_stack_config_id_from_string_yaml={interp_stack_config_id}")
        interp_stack.interp_stack_id = interp_stack_config_id
        return interp_stack

    @classmethod
    def from_string_json(cls, interp_stack_string : str):
        interp_stack = cls.__new__(cls)
        cdef int interp_stack_config_id
        _logger.debug(f"get_interp_stack_config_from_string_json: interp_stack_string={interp_stack_string}")
        yac_cget_interp_stack_config_from_string_json(interp_stack_string.encode(), &interp_stack_config_id)
        _logger.debug(f"interp_stack_config_id_from_string_json={interp_stack_config_id}")
        interp_stack.interp_stack_id = interp_stack_config_id
        return interp_stack

    def add_average(self, reduction_type : AverageReductionType, partial_coverage):
        """
        @see yac_cadd_interp_stack_config_average
        """
        _logger.debug(f"add_interp_stack_config_average interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_average(self.interp_stack_id,
                                             AverageReductionType(reduction_type).value,
                                             partial_coverage)

    def add_ncc(self, reduction_type : NCCReductionType, partial_coverage):
        """
        @see yac_cadd_interp_stack_config_ncc
        """
        _logger.debug(f"add_interp_stack_config_ncc interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_ncc(self.interp_stack_id,
                                         NCCReductionType(reduction_type).value,
                                         partial_coverage)

    def add_nnn(self, reduction_type : NNNReductionType, n : int,
                max_search_distance : _np.float64, scale : _np.float64):
        """
        @see yac_cadd_interp_stack_config_nnn
        """
        _logger.debug(f"add_interp_stack_config_nnn interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_nnn(self.interp_stack_id,
                                         NNNReductionType(reduction_type).value,
                                         n, max_search_distance, scale)

    def add_conservative(self, order : int, enforced_conserv : int,
                         partial_coverage : int, normalisation : ConservNormalizationType):
        """
        @see yac_cadd_interp_stack_config_conservative
        """
        _logger.debug(f"add_interp_stack_config_conservative interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_conservative(self.interp_stack_id,
                                                  order, enforced_conserv,
                                                  partial_coverage,
                                                  ConservNormalizationType(normalisation).value)

    def add_spmap(self, spread_distance : _np.float64, max_search_distance : _np.float64,
                  weight_type : SPMAPWeightType, scale_type : SPMAPScaleType,
                  src_sphere_radius : _np.float64, tgt_sphere_radius : _np.float64):
        """
        @see yac_cadd_interp_stack_config_spmap
        """
        _logger.debug(f"add_interp_stack_config_spmap interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_spmap(self.interp_stack_id,
                                           spread_distance, max_search_distance,
                                           SPMAPWeightType(weight_type).value,
                                           SPMAPScaleType(scale_type).value,
                                           src_sphere_radius, tgt_sphere_radius)

    def add_hcsbb(self):
        """
        @see yac_cadd_interp_stack_config_hcsbb
        """
        _logger.debug(f"add_interp_stack_config_hcsbb interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_hcsbb(self.interp_stack_id)

    def add_user_file(self, filename : str):
        """
        @see yac_cadd_interp_stack_config_user_file
        """
        _logger.debug(f"add_interp_stack_config_user_file interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_user_file(self.interp_stack_id,
                                               filename.encode())

    def add_fixed(self, value : _np.float64):
        """
        @see yac_cadd_interp_stack_config_fixed
        """
        _logger.debug(f"add_interp_stack_config_fixed interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_fixed(self.interp_stack_id, value)

    def add_check(self, constructor_key : str, do_search_key : str):
        """
        @see yac_cadd_interp_stack_config_check
        """
        _logger.debug(f"add_interp_stack_config_check interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_check(self.interp_stack_id,
                                           constructor_key.encode(), do_search_key.encode())

    def add_creep(self, creep_distance : int):
        """
        @see yac_cadd_interp_stack_config_creep
        """
        _logger.debug(f"add_interp_stack_config_creep interp_stack_id={self.interp_stack_id}")
        yac_cadd_interp_stack_config_creep(self.interp_stack_id,
                                           creep_distance)

    def __del__(self):
        """
        @see yac_cfree_interp_stack_config
        """
        _logger.debug(f"free_interp_stack_config interp_stack_id={self.interp_stack_id}")
        yac_cfree_interp_stack_config(self.interp_stack_id)

def version():
    """
    @see yac_cget_version
    """
    return bytes.decode(yac_cget_version())


if Py_AtExit(yac_cfinalize) < 0:
    print(
        b"WARNING: %s\n",
        b"could not register yac_cfinalize with Py_AtExit()",
    )

cdef yac_abort_func _prev_abort_func = yac_get_abort_handler()

cdef void yac_python_abort(MPI_Comm comm, const char* msg,
                           const char* source, int line) noexcept with gil:
    import traceback
    traceback.print_stack()
    _prev_abort_func(comm, msg, source, line)

yac_set_abort_handler(yac_python_abort)
