##
# @file yac.pyx
# @brief cython wrapper for the yac c-API
#
# Contains lightweight cython wrappers for the yac c-API.
#
# @copyright Copyright  (C)  2023 DKRZ, MPI-M
#
# @author Nils-Arne Dreier <dreier@dkrz.de>
#

#
# Keywords:
# Maintainer: Nils-Arne Dreier <dreier@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
#
# This file is part of YAC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import cython
from enum import Enum
import numpy as _np

from libc.stdlib cimport malloc, free

cdef import from "<mpi.h>" nogil:
  ctypedef struct _mpi_comm_t
  ctypedef _mpi_comm_t* MPI_Comm
  int MPI_Comm_c2f(MPI_Comm)
  MPI_Comm MPI_Comm_f2c(int)

cdef extern from "Python.h":
    int Py_AtExit(void (*)())

cdef extern from "yac_interface.h":
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

  cdef const int _NNN_AVG "YAC_NNN_AVG"
  cdef const int _NNN_DIST "YAC_NNN_DIST"
  cdef const int _NNN_GAUSS "YAC_NNN_GAUSS"
  cdef const int _NNN_RBF "YAC_NNN_RBF"

  cdef const int _CONSERV_DESTAREA "YAC_CONSERV_DESTAREA"
  cdef const int _CONSERV_FRACAREA "YAC_CONSERV_FRACAREA"

  cdef const int _SPMAP_AVG "YAC_SPMAP_AVG"
  cdef const int _SPMAP_DIST "YAC_SPMAP_DIST"

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
                               int *grid_id);
  void yac_cdef_points_curve2d ( const int grid_id,
                                 int nbr_points[2],
                                 const int location,
                                 const double *x_points,
                                 const double *y_points,
                                 int *point_id );
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
  void yac_cdef_couple_instance_( int instance_id,
                                  const char * src_comp, const char * src_grid,
                                  const char * src_field, const char * tgt_comp,
                                  const char * tgt_grid, const char * tgt_field,
                                  const char * coupling_timestep, int timeunit,
                                  int time_reduction, int interp_stack_config_id,
                                  int src_lag, int tgt_lag,
                                  const char * weight_file, int mapping_on_source,
                                  double scale_factor, double scale_summand,
                                  int num_src_mask_names,
                                  const char * const * src_mask_names,
                                  const char * tgt_mask_name)
  void yac_cget_ ( const int field_id,
                   const int collection_size,
                   double *recv_field,
                   int    *info,
                   int    *ierror ) nogil
  void yac_cput_ ( const int field_id,
                   const int collection_size,
                   const double * send_field,
                   int *info,
                   int *ierror ) nogil
  void yac_cput_frac_ ( const int field_id,
                        const int collection_size,
                        double *send_field,
                        double *send_frac_mask,
                        int    *info,
                        int    *ierr ) nogil
  const char* yac_cget_field_name_from_field_id ( int field_id )
  int yac_cget_role_from_field_id ( int field_id )
  const char* yac_cget_timestep_from_field_id ( int field_id )
  size_t yac_get_grid_size ( int location, int grid_id )
  size_t yac_get_points_size ( int points_id  )
  int yac_cget_collection_size_from_field_id ( const int field_id )
  const char* yac_cget_field_datetime(int field_id)
  void yac_cget_interp_stack_config(int * interp_stack_config_id)
  void yac_cadd_interp_stack_config_average(
  int interp_stack_config_id, int reduction_type, int partial_coverage)
  void yac_cadd_interp_stack_config_nnn(int interp_stack_config_id, int type,
                                       unsigned int n, double scale)
  void yac_cadd_interp_stack_config_conservative(
      int interp_stack_config_id, int order, int enforced_conserv,
      int partial_coverage, int normalisation)
  void yac_cadd_interp_stack_config_spmap(
      int interp_stack_config_id, double spread_distance,
      double max_search_distance, int weight_type)
  void yac_cadd_interp_stack_config_hcsbb(int interp_stack_config_id)
  void yac_cadd_interp_stack_config_user_file(
      int interp_stack_config_id, char * filename, char * src_grid_name,
      char * tgt_grid_name)
  void yac_cadd_interp_stack_config_fixed(
      int interp_stack_config_id, double value)
  void yac_cadd_interp_stack_config_check(
      int interp_stack_config_id, char * constructor_key, char * do_search_key)
  void yac_cadd_interp_stack_config_creep(
      int interp_stack_config_id, int creep_distance)
  void yac_cfree_interp_stack_config(int interp_stack_config_id)
  void yac_cset_core_mask ( const int * is_core,
                            int location,
                            int grid_id);
  void yac_cset_mask ( const int * is_valid,
                       int points_id )
  void yac_cdef_mask_named ( const int grid_id,
                             const int nbr_points,
                             const int location,
                             const int * is_valid,
                             const char * name,
                             int *mask_id );
  void yac_cfinalize();
  ctypedef void (*yac_abort_func)(MPI_Comm comm, const char *msg,
                                  const char *source, int line) except *
  void yac_set_abort_handler(yac_abort_func custom_abort);
  yac_abort_func yac_get_abort_handler();
  void yac_cread_config_yaml_instance( int yac_instance_id,
                                       const char * yaml_file);

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
                yac_cinit()
                instance_id = yac_cget_default_instance_id()
            else:
                yac_cinit_instance(&instance_id)
        else:
            from mpi4py import MPI
            if type(comm) is MPI.Intracomm:
                comm = MPI.Comm.py2f(comm)
            if default_instance:
                yac_cinit_comm(MPI_Comm_f2c(comm))
                instance_id = yac_cget_default_instance_id()
            else:
                yac_cinit_comm_instance(MPI_Comm_f2c(comm), &instance_id)
        self.instance_id = instance_id
        self.initialized = True

    @classmethod
    @property
    def default_instance(cls):
        yac = cls.__new__(cls)
        yac.instance_id = yac_cget_default_instance_id()
        yac.initialized = False
        return yac

    def __del__(self):
        """
        @see yac_cfinalize_instance
        """
        if self.initialized:
            yac_ccleanup_instance(self.instance_id)


    def def_comp(self, comp_name : str):
        """
        @see yac_cdef_comp_instance
        """
        cdef int comp_id
        yac_cdef_comp_instance(self.instance_id, comp_name.encode(), &comp_id)
        return Component(comp_id)

    def def_comps(self, comp_names = []):
        """
        @see yac_cdef_comps_instance
        """
        cdef int comp_len = len(comp_names)
        cdef const char **c_comp_names = <const char **>malloc(comp_len * sizeof(const char *))
        cdef int *c_comp_ids = <int*>malloc(comp_len * sizeof(int))
        byte_comp_names = [c.encode() for c in comp_names]
        for i in range(comp_len):
            c_comp_names[i] = byte_comp_names[i]
        yac_cdef_comps_instance(self.instance_id, c_comp_names, comp_len, c_comp_ids)
        comp_list = [Component(c_comp_ids[i]) for i in range(comp_len) ]
        free(c_comp_names)
        free(c_comp_ids)
        return comp_list

    def predef_comp(self, comp_name :str):
        """
        @see yac_cpredef_comp_instance
        """
        cdef int comp_id
        yac_cpredef_comp_instance(self.instance_id, comp_name.encode(), &comp_id)
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
            yac_cdef_datetime_instance ( self.instance_id,
                                         start_datetime.encode(),
                                         NULL)

        try:
            if(type(end_datetime) is datetime.datetime):
                end_datetime = end_datetime.isoformat()
        except:
            pass
        if end_datetime is not None:
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
        cdef char ** src_mask_names_ptr = NULL
        cdef char * tgt_mask_name_ptr = NULL
        if src_masks_names is not None:
            if type(src_masks_names) is str:
                src_masks = [src_masks_names]
            src_masks_enc = [s.encode() for s in src_masks_names]
            src_mask_names_ptr = <char**>malloc(len(src_masks_enc) * sizeof(char*))
            for i in range(len(src_masks_enc)):
                src_mask_names_ptr[i] = src_masks_enc[i]
        if tgt_mask_name is not None:
            tgt_mask_enc = tgt_mask_name.encode()
            tgt_mask_name = tgt_mask_enc
        yac_cdef_couple_instance_(self.instance_id,
                                  src_comp.encode(), src_grid.encode(), src_field.encode(),
                                  tgt_comp.encode(), tgt_grid.encode(), tgt_field.encode(),
                                  coupling_timestep.encode(), TimeUnit(timeunit).value,
                                  Reduction(time_reduction).value,
                                  interp_stack.interp_stack_id, src_lag, tgt_lag,
                                  weight_file_ptr, mapping_on_source,
                                  scale_factor, scale_summand,
                                  0, src_mask_names_ptr, tgt_mask_name_ptr)

    def enddef(self):
        """
        @see yac_cenddef_instance
        """
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
        free(ret);
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
        free(ret);
        return grid_list

    def get_comp_grid_names(self, comp_name):
        """
        @see yac_cget_comp_grid_names
        """
        cdef int nbr_grids = yac_cget_comp_nbr_grids_instance(self.instance_id, comp_name.encode())
        cdef const char **ret = <const char **>malloc(nbr_grids * sizeof(const char *))
        yac_cget_comp_grid_names_instance(self.instance_id, comp_name.encode(), nbr_grids, ret)
        grid_list = [bytes(ret[i]).decode('UTF-8') for i in range(nbr_grids) ]
        free(ret);
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
        free(ret);
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
        yac_cenable_field_frac_mask_instance(self.instance_id,
                                             comp_name.encode(),
                                             grid_name.encode(),
                                             field_name.encode(),
                                             frac_mask_fallback_value)

    def def_component_metadata(self, comp_name : str, metadata : bytes):
        """
        @see yac_cset_component_metadata
        """
        yac_cdef_component_metadata_instance(self.instance_id ,
                                             comp_name.encode(), metadata)

    def def_grid_metadata(self, grid_name : str, metadata : bytes):
        """
        @see yac_cset_grid_metadata
        """
        yac_cdef_grid_metadata_instance(self.instance_id,
                                        grid_name.encode(), metadata)

    def def_field_metadata(self, comp_name : str, grid_name : str,
                           field_name : str,metadata : bytes):
        """
        @see yac_cget_field_metadata
        """
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
        yac_cread_config_yaml_instance(self.instance_id, yaml_file.encode())

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
        cdef int[::1] global_index_view = _np.ascontiguousarray(global_index, dtype=_np.intc)
        yac_cset_global_index(&global_index_view[0], location.value, self.grid_id)

    @property
    def nbr_cells(self):
        """
        @see yac_get_grid_size (`int`, read-only)
        """
        return yac_get_grid_size ( Location.CELL.value, self.grid_id )

    @property
    def nbr_corners(self):
        """
        @see yac_get_grid_size (`int`, read-only)
        """
        return yac_get_grid_size ( Location.CORNER.value, self.grid_id )

    @property
    def nbr_edges(self):
        """
        @see yac_get_grid_size (`int`, read-only)
        """
        return yac_get_grid_size ( Location.EDGE.value, self.grid_id )

    def set_core_mask(self, is_core, location : Location):
        """
        @see yac_cset_core_mask
        """
        cdef size_t len_is_core = len(is_core)
        assert(len_is_core == yac_get_grid_size( location.value, self.grid_id ) )
        cdef int[::1] np_mask = _np.ascontiguousarray(is_core, dtype=_np.intc)
        yac_cset_core_mask ( &np_mask[0], location.value, self.grid_id)

    def def_mask(self, location : Location,
                 is_valid, name = None):
        cdef int len_is_valid = len(is_valid)
        cdef int[::1] np_mask = _np.ascontiguousarray(is_valid, dtype=_np.int32)
        cdef int mask_id
        cdef char* c_name = NULL
        if name is not None:
            name_enc = name.encode()
            c_name = name_enc
        yac_cdef_mask_named ( self.grid_id,
                              len_is_valid,
                              location.value,
                              &np_mask[0],
                              c_name,
                              &mask_id );
        return Mask(mask_id)

class Points:
    """
    Stores the points_id and provides further functionality
    """
    def __init__(self, points_id):
        self.points_id = points_id

    @property
    def size(self):
        """
        @see yac_get_points_size (`int`, read-only)
        """
        return yac_get_points_size ( self.points_id )

    def set_mask(self, is_valid):
        """
        @see yac_cset_mask
        """
        cdef size_t len_is_valid = len(is_valid)
        assert len_is_valid==self.size
        cdef int[::1] np_mask = _np.ascontiguousarray(is_valid, dtype=_np.intc)
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
        cdef double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        cdef int[2] cyclic_view = cyclic
        yac_cdef_grid_reg2d(grid_name.encode(), [len(x),len(y)],
                            cyclic_view, &x[0], &y[0], &grid_id)
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_reg2d
        """
        cdef int points_id
        cdef double[::1] x = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef double[::1] y = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        yac_cdef_points_reg2d(self.grid_id, [len(x), len(y)],
                              location.value, &x[0], &y[0], &points_id)
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
        cdef double[::1] x = _np.ascontiguousarray(x_vertices.flatten(), dtype=_np.double)
        cdef double[::1] y = _np.ascontiguousarray(y_vertices.flatten(), dtype=_np.double)
        cdef int[2] cyclic_view = cyclic
        yac_cdef_grid_curve2d(grid_name.encode(),
                              [_np.shape(x_vertices)[1], _np.shape(y_vertices)[0]],
                              cyclic_view, &x[0], &y[0], &grid_id)
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_vertices, y_vertices):
        """
        @see yac_cdef_points_curve2d
        """
        assert x_vertices.shape == y_vertices.shape
        cdef int points_id
        cdef double[::1] x = _np.ascontiguousarray(x_vertices.flatten(), dtype=_np.double)
        cdef double[::1] y = _np.ascontiguousarray(y_vertices.flatten(), dtype=_np.double)
        yac_cdef_points_curve2d(self.grid_id,
                                [_np.shape(x_vertices)[1], _np.shape(x_vertices)[0]],
                                location.value, &x[0], &y[0], &points_id)
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
        cdef int[::1] num_vertices_per_cell_view = _np.ascontiguousarray(num_vertices_per_cell, dtype=_np.intc)
        cdef double[::1] x_vertices_view = _np.ascontiguousarray(x_vertices, dtype=_np.double)
        cdef double[::1] y_vertices_view = _np.ascontiguousarray(y_vertices, dtype=_np.double)
        cdef int[::1] cell_to_vertex_view = _np.ascontiguousarray(cell_to_vertex, dtype=_np.intc)
        if not use_ll_edges:
            yac_cdef_grid_unstruct(grid_name.encode(), len(x_vertices_view),
                                   len(num_vertices_per_cell_view),
                                   &num_vertices_per_cell_view[0],
                                   &x_vertices_view[0], &y_vertices_view[0],
                                   &cell_to_vertex_view[0], &grid_id)
        else:
            yac_cdef_grid_unstruct_ll(grid_name.encode(), len(x_vertices_view),
                                      len(num_vertices_per_cell_view),
                                      &num_vertices_per_cell_view[0],
                                      &x_vertices_view[0], &y_vertices_view[0],
                                      &cell_to_vertex_view[0], &grid_id)
        super().__init__(grid_id)

    def def_points(self, location : Location,
                   x_points, y_points):
        """
        @see yac_cdef_points_unstruct
        """
        cdef int points_id
        cdef double[::1] x_points_view = _np.ascontiguousarray(x_points, dtype=_np.double)
        cdef double[::1] y_points_view = _np.ascontiguousarray(y_points, dtype=_np.double)
        yac_cdef_points_unstruct(self.grid_id, len(x_points_view), location.value,
                                 &x_points_view[0], &y_points_view[0], &points_id)
        return Points(points_id)

class Field:
    """
    Stored the field_id
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
        cdef int[:] point_ids_array = _np.array([p.points_id for p in points], dtype=_np.intc)
        size = points[0].size
        cdef int[:] mask_ids_array
        if masks is None:
            yac_cdef_field(field_name.encode(), comp.comp_id,
                           &point_ids_array[0], len(point_ids_array),
                           collection_size, timestep.encode(), TimeUnit(timeunit).value, &field_id)
        else:
            if not isinstance(masks, Iterable):
                masks = [masks]
            mask_ids_array = _np.array([m.mask_id for m in masks], dtype=_np.intc)
            yac_cdef_field_mask(field_name.encode(), comp.comp_id,
                                &point_ids_array[0], &mask_ids_array[0],
                                len(point_ids_array),
                                collection_size, timestep.encode(), TimeUnit(timeunit).value, &field_id)
        return Field(field_id, size)

    def get(self, buf = None):
        """
        @see yac_cget_

        @param[out] buf    receive buffer, if `None` a numpy array of correct size is allocated
        """
        cdef int info
        cdef int ierror
        if buf is None:
            buf = _np.empty((self.collection_size, self.size), dtype=_np.double)
        buf = _np.ascontiguousarray(buf.reshape(self.collection_size, self.size), dtype=_np.double)
        cdef double[:,::1] buf_view = buf
        cdef int field_id = self.field_id
        cdef int collection_size = self.collection_size
        with cython.nogil:
            yac_cget_(field_id, collection_size, &buf_view[0,0], &info, &ierror)
        if ierror != 0:
            raise RuntimeError("yac_cget returned error number " + str(ierror))
        return buf, Action(info)

    def put(self, buf, frac_mask = None):
        """
        @see yac_cput_
        """
        cdef int info
        cdef int ierror
        cdef double[:,:,::1] buf_view = _np.ascontiguousarray(buf.reshape(-1,self.collection_size,self.size), dtype=_np.double)
        cdef double[:,:,::1] frac_mask_view
        cdef int field_id = self.field_id
        cdef int collection_size = self.collection_size
        if frac_mask is not None:
            frac_mask_view = _np.ascontiguousarray(
                frac_mask.reshape(-1,self.collection_size,self.size), dtype=_np.double)
            with cython.nogil:
                yac_cput_frac_(field_id, collection_size, &buf_view[0,0,0],
                               &frac_mask_view[0,0,0], &info, &ierror)
        else:
            with cython.nogil:
                yac_cput_(field_id, collection_size, &buf_view[0,0,0], &info, &ierror)
        if ierror != 0:
            raise RuntimeError("yac_cput returned error number " + str(ierror))
        return Action(info)

    @property
    def name(self):
        """
        @see yac_cget_field_name_from_field_id
        """
        return bytes.decode(yac_cget_field_name_from_field_id ( self.field_id ))

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

class InterpolationStack:
    def __init__(self):
        """
        @see yac_cget_interp_stack_config
        """
        cdef int interp_stack_config_id
        yac_cget_interp_stack_config(&interp_stack_config_id)
        self.interp_stack_id = interp_stack_config_id

    def add_average(self, reduction_type : AverageReductionType, partial_coverage):
        """
        @see yac_cadd_interp_stack_config_average
        """
        yac_cadd_interp_stack_config_average(self.interp_stack_id,
                                             AverageReductionType(reduction_type).value,
                                             partial_coverage)

    def add_nnn(self, reduction_type : NNNReductionType, n : int, scale : _np.float64):
        """
        @see yac_cadd_interp_stack_config_nnn
        """
        yac_cadd_interp_stack_config_nnn(self.interp_stack_id,
                                         NNNReductionType(reduction_type).value,
                                         n, scale)

    def add_conservative(self, order : int, enforced_conserv : int,
                         partial_coverage : int, normalisation : ConservNormalizationType):
        """
        @see yac_cadd_interp_stack_config_conservative
        """
        yac_cadd_interp_stack_config_conservative(self.interp_stack_id,
                                                  order, enforced_conserv,
                                                  partial_coverage,
                                                  ConservNormalizationType(normalisation).value)

    def add_spmap(self, spread_distance : _np.float64, max_search_distance : _np.float64,
                  weight_type : SPMAPWeightType):
        """
        @see yac_cadd_interp_stack_config_spmap
        """
        yac_cadd_interp_stack_config_spmap(self.interp_stack_id,
                                           spread_distance, max_search_distance,
                                           SPMAPWeightType(weight_type).value)

    def add_hcsbb(self):
        """
        @see yac_cadd_interp_stack_config_hcsbb
        """
        yac_cadd_interp_stack_config_hcsbb(self.interp_stack_id)

    def add_user_file(self, filename : str, src_grid_name : str,
                      tgt_grid_name : str):
        """
        @see yac_cadd_interp_stack_config_user_file
        """
        yac_cadd_interp_stack_config_user_file(self.interp_stack_id,
                                               filename.encode(), src_grid_name.encode(),
                                               tgt_grid_name.encode())

    def add_fixed(self, value : _np.float64):
        """
        @see yac_cadd_interp_stack_config_fixed
        """
        yac_cadd_interp_stack_config_fixed(self.interp_stack_id, value)

    def add_check(self, constructor_key : str, do_search_key : str):
        """
        @see yac_cadd_interp_stack_config_check
        """
        yac_cadd_interp_stack_config_check(self.interp_stack_id,
                                           constructor_key.encode(), do_search_key.encode())

    def add_creep(self, creep_distance : int):
        """
        @see yac_cadd_interp_stack_config_creep
        """
        yac_cadd_interp_stack_config_creep(self.interp_stack_id,
                                           creep_distance)

    def __del__(self):
        """
        @see yac_cfree_interp_stack_config
        """
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
                           const char* source, int line) except *:
    import traceback
    traceback.print_stack()
    _prev_abort_func(comm, msg, source, line)
    #raise Exception(msg)

yac_set_abort_handler(yac_python_abort)
