/**
 * @file yac_interface.h
 * @brief C User Interface
 *
 * An overview of the C and Fortran User Interface can be found
 * \ref interface "here".
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Teresa Holfeld <teresa.holfeld@zmaw.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Teresa Holfeld <teresa.holfeld@zmaw.de>
 * @par Warranty
 *
 * This code has been tested up to a certain level. Defects and weaknesses,
 * which may be included in the code, do not establish any warranties by the
 * authors. The authors do not make any warranty, express or implied, or assume any
 * liability or responsibility for the use, acquisition or application of this
 * software.
 *
 * @par C User API
 *
 * For a detailed description of the C and Fortran interface please
 * see \ref interface.
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Teresa Holfeld <teresa.holfeld@zmaw.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#ifndef YAC_INTERFACE_H
#define YAC_INTERFACE_H

#include <mpi.h>

extern int const YAC_LOCATION_CELL;
extern int const YAC_LOCATION_CORNER;
extern int const YAC_LOCATION_EDGE;

extern int const YAC_EXCHANGE_TYPE_NONE;
extern int const YAC_EXCHANGE_TYPE_SOURCE;
extern int const YAC_EXCHANGE_TYPE_TARGET;

extern int const YAC_ACTION_NONE;
extern int const YAC_ACTION_REDUCTION;
extern int const YAC_ACTION_COUPLING;
extern int const YAC_ACTION_GET_FOR_RESTART;
extern int const YAC_ACTION_PUT_FOR_RESTART;
extern int const YAC_ACTION_OUT_OF_BOUND;

extern int const YAC_REDUCTION_TIME_NONE;
extern int const YAC_REDUCTION_TIME_ACCUMULATE;
extern int const YAC_REDUCTION_TIME_AVERAGE;
extern int const YAC_REDUCTION_TIME_MINIMUM;
extern int const YAC_REDUCTION_TIME_MAXIMUM;

extern int const YAC_TIME_UNIT_MILLISECOND;
extern int const YAC_TIME_UNIT_SECOND;
extern int const YAC_TIME_UNIT_MINUTE;
extern int const YAC_TIME_UNIT_HOUR;
extern int const YAC_TIME_UNIT_DAY;
extern int const YAC_TIME_UNIT_MONTH;
extern int const YAC_TIME_UNIT_YEAR;
extern int const YAC_TIME_UNIT_ISO_FORMAT;

extern int const YAC_CALENDAR_NOT_SET;
extern int const YAC_PROLEPTIC_GREGORIAN;
extern int const YAC_YEAR_OF_365_DAYS;
extern int const YAC_YEAR_OF_360_DAYS;

extern int const YAC_AVG_ARITHMETIC;
extern int const YAC_AVG_DIST;
extern int const YAC_AVG_BARY;

extern int const YAC_NNN_AVG;
extern int const YAC_NNN_DIST;
extern int const YAC_NNN_GAUSS;
extern int const YAC_NNN_RBF;

extern int const YAC_CONSERV_DESTAREA;
extern int const YAC_CONSERV_FRACAREA;

extern int const YAC_SPMAP_AVG;
extern int const YAC_SPMAP_DIST;

extern int YAC_YAML_EMITTER_DEFAULT; //!<  emit to YAML format
extern int YAC_YAML_EMITTER_JSON;    //!<  emit to JSON format

#define YAC_MAX_CHARLEN (132)

/** \example test_dummy_coupling_c.c
 * This example simulates a whole model setup with three components (ocean,
 * atmosphere, io). It uses one process for each component.
 */

/** \example test_dummy_coupling2_c.c
 * This example simulates a whole model setup with two components (ocean,
 * atmosphere) where rank 2 is shared by both components.
 */

/** \example test_dummy_coupling3_c.c
 * This test checks the support for multiple parallel YAC instances.
 * There are three main components and one dummy component. The dummy component
 * does not take part in any coupling. YAC instances are used.
 */

/** \example test_dummy_coupling4_c.c
 * This test checks the various aggregation options.
 */

/** \example test_dummy_coupling5_c.c
 */

/** \example test_dummy_coupling6_c.c
 */

/** \example test_dummy_coupling7_c.c
 * This test checks the frac mask feature.
 */

/** \example test_dummy_coupling8_c.c
 * This test checks the frac mask feature for aggregated fields.
 */

/** \example test_dummy_coupling9_c.c
 * This test checks the setting and usage of named masks.
 */

/** \example test_restart_c.c
 * This example simulates a restart of a
 * coupled configuration with a grid redefiniton.
 */

/** \example test_restart2.c
 * This example simulates a bit more complex restart
 * of a coupled configuration with a grid redefiniton.
 */

/** \example test_mpi_error.c
 * This example tests the yac_mpi_error interface.
 */

/** \example test_multithreading.c
 * This example tests multithreading with YAC.
 */

/** \example test_query_routines_c.c
 * This example tests query interface routines.
 */

/* -----------------------------------------------------------------------------

    C and Fortran description of the user API

   ----------------------------------------------------------------------------- */

/** \example test_mpi_handshake_c.c
 * This examples demonstrates the usage of the MPI handshake algorithm.
 */

/** MPI Handshake (\ref mpi_handshake_detail)

    Splits the provided communicator into group communicators. Each group
    communicator contains all processes of the provided communicator that
    provided the same group name. The order of the group names can be
    arbitrary on each process. A process can be part of multiple groups.

    @param[in]  comm        MPI communicator
    @param[in]  n           number of group communicators to be generated
    @param[in]  group_names group names
    @param[out] group_comms group communicators

    @remark This call is collective for processes in comm.
    @remark \ref yac_cinit will call this routine for MPI_COMM_WORLD
            and provided the group name "yac" to build its internal
            communicator.
    @remark If a process does not wish to use yac but is part of a
            communicator (e.g. MPI_COMM_WORLD) that is used by other
            processes that actually do use yac. This routine can be
            called, while not providing the group name "yac, to exclude
            a process from later collective yac-calls.
*/
void yac_cmpi_handshake(MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms);

/** Getter function for the default instance

    Returns the instance id of the default instance.

    @remark The default instance needs to be initialised before calling this function
    @returns The instance id of the default isntance

 */
int yac_cget_default_instance_id();

/** Elementary initialisation of the whole system

     @remark this call initialises the default YAC instance
     @remark This call executes a MPI handshake (\ref mpi_handshake_detail) on
     MPI_COMM_WORLD with the group name "yac" and uses the created communicator
     for the initialization. If you dont want to execute a MPI Handshake on
     MPI_COMM_WORLD you can call yac_cinit_comm with MPI_COMM_WORLD directly.
*/
    void yac_cinit (void);

/** Elementary initialisation of the whole system

     @param[out] yac_instance_id id of the YAC instance initialised by
                                 this call
     @remark This call executes a MPI handshake (\ref mpi_handshake_detail) on
     MPI_COMM_WORLD with the group name "yac" and uses the created communicator
     for the initialization. If you dont want to execute a MPI Handshake on
     MPI_COMM_WORLD you can call yac_cinit_comm_instance with MPI_COMM_WORLD directly.
*/
    void yac_cinit_instance ( int * yac_instance_id );

/** Elementary initialisation of the whole system using a user-provided
 *  MPI communicator

     @param[in] comm            MPI communicator
     @remark this call initialises the default YAC instance
     @remark the MPI communicator provided to this routine has to contain
             all processes that will take part in the coupling
     @remark this call is collective for all processes in comm
*/
    void yac_cinit_comm ( MPI_Comm comm );

/** Elementary initialisation of the whole system using a user-provided
 *  MPI communicator

     @param[in]  comm            MPI communicator
     @param[out] yac_instance_id id of the YAC instance initialised by
                                 this call
     @remark In case the user has multiple YAC instances in parallel, he has
             to initialise YAXT himself.
*/
    void yac_cinit_comm_instance ( MPI_Comm comm, int * yac_instance_id );

/* -------------------------------------------------------------------------------- */

/** Dummy for initialisation of the whole system
     @remark This routine can be called instead of \ref yac_cinit or
             \ref yac_cinit_instance, if the local process does not wish to
             initialise a YAC instance.
*/
    void yac_cinit_dummy(void);

/** Dummy initialisation of the whole system using a user-provided
 *  MPI world communicator

     @param[in] comm MPI communicator
     @remark This routine can be called instead of \ref yac_cinit_comm or
             \ref yac_cinit_comm_instance, if the local process does not wish to
             initialise a YAC instance.
*/
    void yac_cinit_comm_dummy ( MPI_Comm comm );

/* -------------------------------------------------------------------------------- */

/** Function for reading configuration from YAML configuration file for a specific
 *  YAC instance

     @param[in] yac_instance_id YAC instance id
     @param[in] yaml_file       YAML configuration file
*/
void yac_cread_config_yaml_instance( int yac_instance_id,
  const char * yaml_file);

/** Function for reading configuration from YAML configuration file

     @param[in] yaml_file       YAML configuration file
*/
void yac_cread_config_yaml(const char * yaml_file);

/** Function for reading configuration from JSON configuration file for a specific
 *  YAC instance

     @param[in] yac_instance_id YAC instance id
     @param[in] json_file       JSON configuration file
*/
void yac_cread_config_json_instance( int yac_instance_id,
  const char * json_file);

/** Function for reading configuration from JSON configuration file

     @param[in] json_file       JSON configuration file
*/
void yac_cread_config_json(const char * json_file);

/* -------------------------------------------------------------------------------- */

/** Clean-up YAC (same as yac_cfinalize but without MPI_Finalize; see
 *  \ref phase_restart)
     @remark In case the user initialised yaxt himself, he has to finalise it after
             yac_ccleanup in order to be able to call yac_cinit again.
     @remark cleans up default YAC instance
*/
     void yac_ccleanup ();

/** Clean-up YAC (same as yac_cfinalize_instance but without MPI_Finalize;
 *  see \ref phase_restart)
     @param[in] yac_instance_id id of the YAC instance to be cleaned up
     @remark In case the user initialised yaxt himself, he has to finalise it after
             yac_ccleanup in order to be able to call yac_cinit again.
*/
     void yac_ccleanup_instance ( int yac_instance_id );

/* -------------------------------------------------------------------------------- */

/** Cleans up default YAC instance (if available) and finales MPI,
 *  if it was initialised by YAC
*/
     void yac_cfinalize ();

/** Cleans up provided YAC instance and finales MPI,
 *  if it was initialised by YAC
     @param[in] yac_instance_id id of the YAC instance to be finalised
*/
     void yac_cfinalize_instance (int yac_instance_id);

/* -------------------------------------------------------------------------------- */

/** Definition of job start and end datetime for the default YAC instance

     @param[in] start_datetime calendar job start datetime
     @param[in] end_datetime   calendar job end datetime
 */

     void yac_cdef_datetime ( const char * start_datetime,
                              const char * end_datetime );

/** Definition of job start and end datetime

     @param[in] yac_instance_id id of the YAC instance
     @param[in] start_datetime  calendar job start datetime
     @param[in] end_datetime    calendar job end datetime
 */

     void yac_cdef_datetime_instance ( int yac_instance_id,
                                       const char * start_datetime,
                                       const char * end_datetime );

/** Defines the calendar of the default instance

    @param[in] calendar The calendar @see YAC_PROLEPTIC_GREGORIAN etc.
*/

     void yac_cdef_calendar ( int calendar );

/* -------------------------------------------------------------------------------- */

/** Providing a group communicator for the default YAC instance

     @param[out] group_comm communicator containing all processes that
                            passed the same group name to the init-routine
                            or MPI_COMM_NULL if not group name was provided.
 */

     void yac_cget_groupcomm( MPI_Comm * group_comm );

/** Providing a group communicator

     @param[in] yac_instance_id id of the YAC instance
     @param[out] group_comm communicator containing all processes that
                            passed the same group name to the init-routine
                            or MPI_COMM_NULL if not group name was provided.
 */

     void yac_cget_groupcomm_instance( int yac_instance_id,
                                       MPI_Comm * group_comm );

/* -------------------------------------------------------------------------------- */

/** Non-blocking pre-definition of a component for the default YAC instance.

     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call does not finish the component defintiion phase. I.e. the
     component communicators are not set up after this call.
*/

void yac_cpredef_comp(char const * comp_name, int * comp_id);

/** Non-blocking pre-definition of a component.

    @param[in]  yac_instance_id id of the YAC instance
    @param[in]  comp_name name of the component
    @param[out] comp_id   component Id
    @remark this call does not finish the component defintiion phase. I.e. the
    component communicators are not set up after this call.
*/

void yac_cpredef_comp_instance(int yac_instance_id,
                               char const * comp_name,
                               int * comp_id);

/* -------------------------------------------------------------------------------- */

/** Elementary definition of the component for the default YAC instance

     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call is collective for all processes that initialised the
             default YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance. If you want to predefine a component, please use `yac_cpredef_comp`.
*/

     void yac_cdef_comp ( const char * comp_name,
                          int * comp_id );

/** Elementary definition of the component.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call is collective for all processes that initialised the
             provided YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance. If you want to predefine a component, please use
             `yac_cpredef_comp_instance`.
*/

     void yac_cdef_comp_instance ( int yac_instance_id,
                                   const char * comp_name,
                                   int * comp_id );

/* -------------------------------------------------------------------------------- */

/** Elementary definition of the components for the default YAC instance

     @param[in]  comp_names names of the components
     @param[in]  num_comps  number of components
     @param[out] comp_ids   component Id's
     @remark this call is collective for all processes that initialised the
             default YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance
*/

     void yac_cdef_comps ( const char ** comp_names,
                           int num_comps,
                           int * comp_ids );

/** Elementary definition of the components.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_names names of the components
     @param[in]  num_comps  number of components
     @param[out] comp_ids   component Id's
     @remark this call is collective for all processes that initialised the
             provided YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance
*/

     void yac_cdef_comps_instance ( int yac_instance_id,
                                    const char ** comp_names,
                                    int num_comps,
                                    int * comp_ids );

/* -------------------------------------------------------------------------------- */

/** Providing the component MPI communicator

     @param[in]  comp_id   local component ID (from yac_cdef_comp),
                           handle to the component struct
     @param[out] comp_comm component communicator

*/
     void yac_cget_comp_comm ( int comp_id,
                               MPI_Comm *comp_comm );

/* -------------------------------------------------------------------------------- */

/** Generates an MPI communicator that contains all processes of the provided
 *  components

     @param[in]  comp_names name of components
     @param[in]  num_comps  number of components
     @param[out] comps_comm generated communicator
     @remark the components have to defined in the default YAC instance
     @remark the local process has to be in at least one of the provided
             components
     @remark this call is collective for all processes in the provided
             components
*/
     void yac_cget_comps_comm ( const char ** comp_names,
                                int num_comps,
                                MPI_Comm * comps_comm);

/** Generates an MPI communicator that contains all processes of the provided
 *  components

     @param[in]  yac_instance_id    id of the YAC instance
     @param[in]  comp_names name of components
     @param[in]  num_comps  number of components
     @param[out] comps_comm generated communicator
     @remark the local process has to be in at least one of the provided
             components
     @remark this call is collective for all processes in the provided
             components
*/
     void yac_cget_comps_comm_instance ( int yac_instance_id,
                                         const char ** comp_names,
                                         int num_comps,
                                         MPI_Comm * comps_comm );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D regular grids.

     @param[in]  grid_id    handle to the grid struct (returned from yac_cdef_grid)
     @param[in]  nbr_points number of points
     @param[in]  location   location of points [edge|vertex|center]
     @param[in]  x_points   array of point longitudes (in rad)
     @param[in]  y_points   array of point latitudes (in rad)
     @param[out] point_id   returned handle to the points struct
*/

     void yac_cdef_points_reg2d ( int const grid_id,
                                  int const *nbr_points,
                                  int const location,
                                  double const *x_points,
                                  double const *y_points,
                                  int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D curvilinear grids.

     @param[in]  grid_id    handle to the grid struct (returned from yac_cdef_grid)
     @param[in]  nbr_points number of points
     @param[in]  location   location of points [edge|vertex|center]
     @param[in]  x_points   array of point longitudes (in rad)
     @param[in]  y_points   array of point latitudes (in rad)
     @param[out] point_id   returned handle to the points struct
     @remark The array x_points and y_points are expected to be of size
             nbr_points[0] * nbr_points[1].
*/

     void yac_cdef_points_curve2d ( int const grid_id,
                                    int const *nbr_points,
                                    int const location,
                                    double const *x_points,
                                    double const *y_points,
                                    int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D unstructured grids.

     @param[in]  grid_id    handle to the grid struct (returned from yac_cdef_grid)
     @param[in]  nbr_points number of points
     @param[in]  location   location of points [edge|vertex|center]
     @param[in]  x_points   array of point longitudes (in rad)
     @param[in]  y_points   array of point latitudes (in rad)
     @param[out] point_id   returned handle to the points struct
*/

     void yac_cdef_points_unstruct ( int const grid_id,
                                     int const nbr_points,
                                     int const location,
                                     double const *x_points,
                                     double const *y_points,
                                     int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a 2d regular grid

     @param[in]  grid_name    name of the grid
     @param[in]  nbr_vertices 2d array containing the number of vertices in each dimension
     @param[in]  cyclic       2d array containing information about cyclic behaviour in each dimension
     @param[in]  x_vertices   array of vertex longitudes (in rad)
     @param[in]  y_vertices   array of vertex latitudes (in rad)
     @param[out] grid_id      id of generated grid
     @remark on the local order of edges:
             generate for each edge a tuple containing the local indices of the
             adjacent vertices (lowest index first) and sort this list by the
             vertex indices; this gives the the local order of edges
*/
     void yac_cdef_grid_reg2d ( const char * grid_name,
                                int nbr_vertices[2],
                                int cyclic[2],
                                double *x_vertices,
                                double *y_vertices,
                                int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of a 2d curvilinear grid

     @param[in]  grid_name    name of the grid
     @param[in]  nbr_vertices 2d array containing the number of vertices in each dimension
     @param[in]  cyclic       2d array containing information about cyclic behaviour in each dimension
     @param[in]  x_vertices   array of vertex longitudes (in rad)
     @param[in]  y_vertices   array of vertex latitudes (in rad)
     @param[out] grid_id      id of generated grid
     @remark on the local order of edges:
             generate for each edge a tuple containing the local indices of the
             adjacent vertices (lowest index first) and sort this list by the
             vertex indices; this gives the the local order of edges
     @remark The array x_vertices and y_vertices are expected to be of size
             nbr_vertices[0] * nbr_vertices[1].
*/
     void yac_cdef_grid_curve2d ( const char * grid_name,
                                  int nbr_vertices[2],
                                  int cyclic[2],
                                  double *x_vertices,
                                  double *y_vertices,
                                  int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid

     @param[in]  grid_name             name of the grid
     @param[in]  nbr_vertices          number of vertices
     @param[in]  nbr_cells             number of cells
     @param[in]  num_vertices_per_cell array containing the number of vertices for each cell
     @param[in]  x_vertices            array of vertex longitudes (in rad)
     @param[in]  y_vertices            array of vertex latitudes (in rad)
     @param[in]  cell_to_vertex        connectivity of vertices belonging to cells\n
                                       (the vertex indices per cell have to be in clockwise or counterclockwise ordering)
     @param[out] grid_id               id of generated grid
*/
     void yac_cdef_grid_unstruct ( const char * grid_name,
                                   int nbr_vertices,
                                   int nbr_cells,
                                   int *num_vertices_per_cell,
                                   double *x_vertices,
                                   double *y_vertices,
                                   int *cell_to_vertex,
                                   int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid with lon-lat edges

     @param[in]  grid_name             name of the grid
     @param[in]  nbr_vertices          number of vertices
     @param[in]  nbr_cells             number of cells
     @param[in]  num_vertices_per_cell array containing the number of vertices for each cell
     @param[in]  x_vertices            array of vertex longitudes (in rad)
     @param[in]  y_vertices            array of vertex latitudes (in rad)
     @param[in]  cell_to_vertex        connectivity of vertices belonging to cells\n
                                       (the vertex indices per cell have to be in clockwise or counterclockwise ordering)
     @param[out] grid_id               id of generated grid
     @remark YAC will check all edges of the grid and determine whether they are on
             circles of longitudes (same x coordinate) or latitudes (same y coordinate).
             An edge that does not fullfill this condition will cause an error.
*/
     void yac_cdef_grid_unstruct_ll ( const char * grid_name,
                                      int nbr_vertices,
                                      int nbr_cells,
                                      int *num_vertices_per_cell,
                                      double *x_vertices,
                                      double *y_vertices,
                                      int *cell_to_vertex,
                                      int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Set global ids for a grid

     @param[in] global_index array of global indices
     @param[in] location     cell/vertex/edge
     @param[in] grid_id      grid id
     @remark global indices are to be provided in the local order of their
             respective cells/vertices/edges
*/

     void yac_cset_global_index ( int const * global_index,
                                  int location,
                                  int grid_id);

/* -------------------------------------------------------------------------------- */

/** Set core mask for a grid

     @param[in] is_core   0 for cells/vertices/edges that are halos,
                          1 for cells/vertices/edges that are core
     @param[in] location  cell/vertex/edge
     @param[in] grid_id   grid id
     @remark cells/vertices/edges who are halos are not used or set
             in a put or get operation
     @remark core mask values are to be provided in the local order of their
             respective cells/vertices/edges
*/

     void yac_cset_core_mask ( int const * is_core,
                               int location,
                               int grid_id);

/* -------------------------------------------------------------------------------- */

/** Set the default mask for points

     @param[in] is_valid   0 for points that are masked out, 1 for valid points
     @param[in] points_id  id of points/cells
*/

     void yac_cset_mask ( int const * is_valid,
                          int points_id );

/* -------------------------------------------------------------------------------- */

/** define a mask for a grid

     @param[in]  grid_id    grid ID
     @param[in]  nbr_points number of points
     @param[in]  location   cell/vertex/edge
     @param[in]  is_valid   0 for points that are masked out, 1 for valid points
     @param[out] mask_id    returned handle to the mask struct
*/

     void yac_cdef_mask ( int const grid_id,
                          int const nbr_points,
                          int const location,
                          int const * is_valid,
                          int *mask_id );

/* -------------------------------------------------------------------------------- */

/** define a named mask for a grid

     @param[in]  grid_id    grid ID
     @param[in]  nbr_points number of points
     @param[in]  location   cell/vertex/edge
     @param[in]  is_valid   0 for points that are masked out, 1 for valid points
     @param[in]  name       name of the mask
     @param[out] mask_id    returned handle to the mask struct
*/

     void yac_cdef_mask_named ( int const grid_id,
                                int const nbr_points,
                                int const location,
                                int const * is_valid,
                                char const * name,
                                int *mask_id );

/* -------------------------------------------------------------------------------- */

/** Definition of the coupling field (using default mask, if defined)

     @param[in]  field_name      character string providing the name of the coupling field
     @param[in]  component_id    component ID
     @param[in]  point_ids       point IDs
     @param[in]  num_pointsets   number of pointsets per grid
     @param[in]  collection_size collection size
     @param[in]  timestep        timestep
     @param[in]  time_unit       time unit
     @param[out] field_id        returned field_id which has to be used to identify coupling fields
                                 in yac_cput and yac_cget.
     @remark A calendar has to be defined before this routine is called.
             This can be done by a call to \ref yac_cdef_calendar, by reading
             of configuration file that defines a calendar, or by
             synchronizing with other processes that already have defined
             a calendar using \ref yac_csync_def.
*/

       void yac_cdef_field ( char const * field_name,
                             int const component_id,
                             int const * point_ids,
                             int const num_pointsets,
                             int collection_size,
                             const char* timestep,
                             int time_unit,
                             int * field_id );

/* -------------------------------------------------------------------------------- */

/** Definition of the coupling field

     @param[in]  field_name      character string providing the name of the coupling field
     @param[in]  component_id    component ID
     @param[in]  point_ids       point IDs
     @param[in]  mask_ids        mask IDs
     @param[in]  num_pointsets   number of pointsets per grid
     @param[in]  collection_size collection size
     @param[in]  timestep        timestep
     @param[in]  time_unit
     @param[out] field_id        returned field_id which has to be used to identify coupling fields
                                 in yac_cput and yac_cget.
     @remark A calendar has to be defined before this routine is called.
             This can be done by a call to \ref yac_cdef_calendar, by reading
             of configuration file that defines a calendar, or by
             synchronizing with other processes that already have defined
             a calendar using \ref yac_csync_def.
*/

       void yac_cdef_field_mask ( char const * field_name,
                                  int const component_id,
                                  int const * point_ids,
                                  int const * mask_ids,
                                  int const num_pointsets,
                                  int collection_size,
                                  const char* timestep,
                                  int time_unit,
                                  int * field_id );

/* -------------------------------------------------------------------------- */

/** Enables fractional masking for a coupling field and sets the fallback value

     @param[in]  comp_name
     @param[in]  grid_name
     @param[in]  field_name
     @param[in]  frac_mask_fallback_value fractional mask fallback value
     @remark If this field is used as a source in couple, the user will have
             to provide a fractional mask for all source points along with the
             source field values in all put/exchange operations of this
             field. (see \ref frac_mask_desc)
*/

     void yac_cenable_field_frac_mask(
       const char* comp_name, const char* grid_name, const char* field_name,
       double frac_mask_fallback_value);

/** Enables fractional masking for a coupling field and sets the fallback value

    @param[in] yac_instance_id
    @param[in]  comp_name
    @param[in]  grid_name
    @param[in]  field_name
    @param[in]  frac_mask_fallback_value fractional mask fallback value
    @remark If this field is used as a source in couple, the user will have
            to provide a fractional mask for all source points along with the
            source field values in all put/exchange operations of this
            field. (see \ref frac_mask_desc)
*/

     void yac_cenable_field_frac_mask_instance(
       int yac_instance_id, const char* comp_name, const char* grid_name,
       const char* field_name, double frac_mask_fallback_value);

/* -------------------------------------------------------------------------- */

/** Define metadata for a component
    @param[in] comp_name
    @param[in] metadata
*/
void yac_cdef_component_metadata(const char* comp_name, const char* metadata);

/** Define metadata for a component
    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] metadata
*/
void yac_cdef_component_metadata_instance(int yac_instance_id, const char* comp_name,
  const char* metadata);

/** Define metadata for a grid for the default YAC instance
    @param[in] grid_name
    @param[in] metadata
*/
void yac_cdef_grid_metadata(const char* grid_name, const char* metadata);

/** Define metadata for a grid
    @param[in] yac_instance_id
    @param[in] grid_name
    @param[in] metadata
*/
void yac_cdef_grid_metadata_instance(int yac_instance_id, const char* grid_name,
  const char* metadata);

/** Define metadata for a field
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @param[in] metadata
*/
void yac_cdef_field_metadata(const char* comp_name, const char* grid_name,
  const char* field_name, const char* metadata);

/** Define metadata for a field
    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @param[in] metadata
*/
void yac_cdef_field_metadata_instance(int yac_instance_id,const char* comp_name,
  const char* grid_name, const char* field_name, const char* metadata);

/* -------------------------------------------------------------------------- */

/** Returns the action a put/get would return for the current timestep

@param[in]  field_id field id returned by \ref yac_cdef_field
@param[out] action   action for the current timestep
                          (\ref YAC_ACTION_NONE,
                           \ref YAC_ACTION_REDUCTION,
                           \ref YAC_ACTION_COUPLING,
                           \ref YAC_ACTION_GET_FOR_RESTART,
                           \ref YAC_ACTION_PUT_FOR_RESTART, or
                           \ref YAC_ACTION_OUT_OF_BOUND)
@remark If the exchange type of the field is \ref YAC_EXCHANGE_TYPE_NONE,
        this routine will return \ref YAC_ACTION_NONE.
@remark If this routine is called before \ref yac_cenddef, the behaviour
        is undefined.
*/

void yac_cget_action(int field_id, int * action);

/** Returns the current datetime of the field
 * @param[in] field_id
 * @returns current datetime of the field
 */
const char * yac_cget_field_datetime(int field_id);

/* -------------------------------------------------------------------------- */

/** If the action for the current timestep is \ref YAC_ACTION_NONE,
    this routine will advance the internal clock to next timestep.

@param[in]  field_id field id returned by \ref yac_cdef_field
@remark If this routine is called before \ref yac_cenddef, the behaviour
        is undefined.
@remark This routine can be called instead of \ref yac_cput or \ref yac_cget,
        if these routines would not have taken any action for the current
        timestep, in order to advance the internal clock.
*/

void yac_cupdate(int field_id);

/* -------------------------------------------------------------------------- */

/** Get an extended coupling configuration (all parameters are set
    to the default values)
    @param[out] ext_couple_config_id extended coupling configuration
 */
void yac_cget_ext_couple_config(int * ext_couple_config_id);

/** Frees an extended coupling configuration
    @param[in] ext_couple_config_id
 */
void yac_cfree_ext_couple_config(int ext_couple_config_id);

/** Sets the weight file name
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] weight_file name of a weight file
    @remark if (weight_file == NULL) no weight file will be written
    @remark YAC will write out the wight file in parallel.
            This parallel output can be configured as described here:
            \ref io_config_detail
 */
void yac_cset_ext_couple_config_weight_file(int ext_couple_config_id,
                                            char const * weight_file);

/** Gets the weight file name
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] weight_file          name of the weight file, if it is set;
                                     NULL otherwise
 */
void yac_cget_ext_couple_config_weight_file(int ext_couple_config_id,
                                            char const ** weight_file);

/** Sets the mapping side
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] mapping_side mapping side
    @remark Valid values are: source = 1 and target = 0
 */
void yac_cset_ext_couple_config_mapping_side(int ext_couple_config_id,
                                             int mapping_side);

/** Gets the mapping side
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] mapping_side         mapping side
 */
void yac_cget_ext_couple_config_mapping_side(int ext_couple_config_id,
                                             int * mapping_side);

/** Sets the scale factor
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] scale_factor         scale factor
 */
void yac_cset_ext_couple_config_scale_factor(int ext_couple_config_id,
                                             double scale_factor);

/** Gets the scale factor
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] scale_factor         scale_factor
 */
void yac_cget_ext_couple_config_scale_factor(int ext_couple_config_id,
                                             double * scale_factor);

/** Sets the scale summand
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] scale_summand        scale summand
 */
void yac_cset_ext_couple_config_scale_summand(int ext_couple_config_id,
                                              double scale_summand);

/** Gets the scale summand
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] scale_summand        scale_summand
 */
void yac_cget_ext_couple_config_scale_summand(int ext_couple_config_id,
                                              double * scale_summand);

/** Sets source mask names
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] num_src_mask_names   number of source mask names
    @param[in] src_mask_names       source mask names
 */
void yac_cset_ext_couple_config_src_mask_names(
  int ext_couple_config_id, size_t num_src_mask_names,
  char const * const * src_mask_names);

/** Gets source mask names
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] num_src_mask_names   number of source mask names
    @param[out] src_mask_names       source mask names
 */
void yac_cget_ext_couple_config_src_mask_names(
  int ext_couple_config_id, size_t * num_src_mask_names,
  char const * const ** src_mask_names);

/** Sets target mask name
    @param[in] ext_couple_config_id extended coupling configuration
    @param[in] tgt_mask_name        target mask name
 */
void yac_cset_ext_couple_config_tgt_mask_name(
  int ext_couple_config_id, char const * tgt_mask_name);

/** Gets target mask name
    @param[in]  ext_couple_config_id extended coupling configuration
    @param[out] tgt_mask_name        target mask name
 */
void yac_cget_ext_couple_config_mask_name(
  int ext_couple_config_id, char const ** tgt_mask_name);

/* -------------------------------------------------------------------------- */

/** Define couple using an extended coupling configuration
    @param[in] yac_instance_id id of the YAC instance for that the couple
               should be defined
    @param[in] src_comp_name component name of the source component
    @param[in] src_grid_name grid name of the source grid
    @param[in] src_field_name field name of the source field
    @param[in] tgt_comp_name component name of the target component
    @param[in] tgt_grid_name grid name of the target grid
    @param[in] tgt_field_name field name of the target field
    @param[in] coupling_timestep time step for the coupling
    @param[in] time_unit unit of coupling_timestep argument
    @param[in] time_reduction type for reducing multiple timesteps
               @see YAC_REDUCTION_TIME_NONE etc.
    @param[in] interp_stack_config_id id of the interpolation stack config to
               be used
    @param[in] src_lag lag for this couple on the source component
    @param[in] tgt_lag lag for this couple on the target component
    @param[in] ext_couple_config_id extended coupling configuration
 */

void yac_cdef_couple_custom_instance(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  int ext_couple_config_id);

/** Define couple in the default yac instance using an extended
    coupling configuration
    @param[in] src_comp_name component name of the source component
    @param[in] src_grid_name grid name of the source grid
    @param[in] src_field_name field name of the source field
    @param[in] tgt_comp_name component name of the target component
    @param[in] tgt_grid_name grid name of the target grid
    @param[in] tgt_field_name field name of the target field
    @param[in] coupling_timestep time step for the coupling given in the defined
               timestep unit
    @param[in] time_unit unit of coupling_timestep argument
    @param[in] time_reduction type for reducing multiple timesteps
               @see YAC_REDUCTION_TIME_NONE etc.
    @param[in] interp_stack_config_id id of the interpolation stack config to
               be used
    @param[in] src_lag lag for this couple on the source component
    @param[in] tgt_lag lag for this couple on the target component
    @param[in] ext_couple_config_id extended coupling configuration
 */

void yac_cdef_couple_custom(
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  int ext_couple_config_id);

/* -------------------------------------------------------------------------- */

/** Define couple
    @param[in] yac_instance_id id of the YAC instance for that the couple
               should be defined
    @param[in] src_comp_name component name of the source component
    @param[in] src_grid_name grid name of the source grid
    @param[in] src_field_name field name of the source field
    @param[in] tgt_comp_name component name of the target component
    @param[in] tgt_grid_name grid name of the target grid
    @param[in] tgt_field_name field name of the target field
    @param[in] coupling_timestep time step for the coupling
    @param[in] time_unit unit of coupling_timestep argument
    @param[in] time_reduction type for reducing multiple timesteps
               @see YAC_REDUCTION_TIME_NONE etc.
    @param[in] interp_stack_config_id id of the interpolation stack config to
               be used
    @param[in] src_lag lag for this couple on the source component
    @param[in] tgt_lag lag for this couple on the target component
 */

void yac_cdef_couple_instance(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag);

/** Define couple in the default yac instance
    @param[in] src_comp_name component name of the source component
    @param[in] src_grid_name grid name of the source grid
    @param[in] src_field_name field name of the source field
    @param[in] tgt_comp_name component name of the target component
    @param[in] tgt_grid_name grid name of the target grid
    @param[in] tgt_field_name field name of the target field
    @param[in] coupling_timestep time step for the coupling given in the defined
               timestep unit
    @param[in] time_unit unit of coupling_timestep argument
    @param[in] time_reduction type for reducing multiple timesteps
               @see YAC_REDUCTION_TIME_NONE etc.
    @param[in] interp_stack_config_id id of the interpolation stack config to
               be used
    @param[in] src_lag lag for this couple on the source component
    @param[in] tgt_lag lag for this couple on the target component
 */

void yac_cdef_couple(
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag);


/* -------------------------------------------------------------------------------- */

/** Checks whether the field dimensions match with field definition and aborts if
 *  there is a mismatch

     @param[in] field_id
     @param[in] collection_size collection size
     @param[in] num_pointsets   number of pointsets
     @param[in] pointset_sizes  size of each pointset
     @remark if num_pointsets is -1, its value will not be checked
     @remark if pointset_sizes is NULL, its values will not be checked
*/
     void yac_ccheck_field_dimensions ( int field_id,
                                        int collection_size,
                                        int num_pointsets,
                                        int const * pointset_sizes );

/* -------------------------------------------------------------------------------- */

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[in]  info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
*/
     void yac_cget_ ( int const field_id,
                      int const collection_size,
                      double *recv_field,
                      int    *info,
                      int    *ierror );

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
*/
     void yac_cget ( int const field_id,
                     int const collection_size,
                     double **recv_field,
                     int    *info,
                     int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_ ( int const field_id,
                      int const collection_size,
                      double *send_field,
                      int    *info,
                      int    *ierror );

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput ( int const field_id,
                     int const collection_size,
                     double *** const send_field,
                     int *info,
                     int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask    - send fractional mask (all data is stored in one
                                                           contiguous part of the memory)
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac_ ( int const field_id,
                           int const collection_size,
                           double *send_field,
                           double *send_frac_mask,
                           int    *info,
                           int    *ierror );

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac ( int const field_id,
                          int const collection_size,
                          double *** const send_field,
                          double *** const send_frac_mask,
                          int *info,
                          int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_ptr_ ( int const field_id,
                          int const collection_size,
                          double ** send_field,
                          int    *info,
                          int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size * nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac_ptr_ ( int const field_id,
                               int const collection_size,
                               double ** send_field,
                               double ** send_frac_mask,
                               int    *info,
                               int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_ ( int const send_field_id,
                           int const recv_field_id,
                           int const collection_size,
                           double *send_field ,
                           double *recv_field,
                           int    *send_info,
                           int    *recv_info,
                           int    *ierror );


/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange ( int const send_field_id,
                          int const recv_field_id,
                          int const collection_size,
                          double *** const send_field,
                          double ** recv_field,
                          int *send_info,
                          int *recv_info,
                          int *ierror );

/* -------------------------------------------------------------------------------- */

/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_frac_ ( int const send_field_id,
                                int const recv_field_id,
                                int const collection_size,
                                double *send_field,
                                double *send_frac_mask,
                                double *recv_field,
                                int    *send_info,
                                int    *recv_info,
                                int    *ierror );


/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_frac ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double *** const send_field,
                               double *** const send_frac_mask,
                               double ** recv_field,
                               int *send_info,
                               int *recv_info,
                               int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_ptr_ ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double ** send_field,
                               double ** recv_field,
                               int    *send_info,
                               int    *recv_info,
                               int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send buffer
                                      dimensions send_frac_mask[collection_size * nbr_fields]
                                                               [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_frac_ptr_ ( int const send_field_id,
                                    int const recv_field_id,
                                    int const collection_size,
                                    double ** send_field,
                                    double ** send_frac_mask,
                                    double ** recv_field,
                                    int    *send_info,
                                    int    *recv_info,
                                    int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Determines whether there is an asynchronous communication associated with
    a field, which is not yet completed (for example by a previous put)

     @param[in]   field_id
     @param[out]  flag     "0" if there is an uncompleted asynchronous
                           communication associated with the field,
                           "1" otherwise
*/

     void yac_ctest ( int field_id, int * flag );

/* -------------------------------------------------------------------------------- */

/** synchronize grids and fields of the default instance
 */
void yac_csync_def ( void );

/** synchronize grids and fields
 *
 *  @param[in]  yac_instance_id id of the YAC instance
 */
void yac_csync_def_instance ( int yac_instance_id );

/** End of the definition phase for the default YAC instance and
 *  setup of internal data structers for coupling (communication matricies
 *  and interpolation weights)
 */
void yac_cenddef ( void );

/** End of the definition phase, invocation of the search and
 *  setup of internal data structers for coupling (communication matricies
 *  and interpolation weights)
 *
 *  @param[in]  yac_instance_id id of the YAC instance
 */
void yac_cenddef_instance ( int yac_instance_id );

/** End of the definition phase for the default YAC instance,
 *  invocation of the search and
 *  setup of internal data structers for coupling (communication matricies
 *  and interpolation weights). In addition, the collected coupling configuration
 *  from all processes is emitted.
 *
 *  @param[in]  emit_flags flags for configuring the generated coupling
 *                         configuration output
 *                         (\ref YAC_YAML_EMITTER_DEFAULT or
 *                          \ref YAC_YAML_EMITTER_JSON)
 *  @param[out] config     coupling configuration
 *  @remark the user is responsible for freeing the returned coupling
 *          configuration
 */

void yac_cenddef_and_emit_config ( int emit_flags, char ** config );

/** End of the definition phase, invocation of the search and
 *  setup of internal data structers for coupling (communication matricies
 *  and interpolation weights). In addition, the collected coupling configuration
 *  from all processes is emitted.
 *
 *  @param[in]  yac_instance_id id of the YAC instance
 *  @param[in]  emit_flags flags for configuring the generated coupling
 *                         configuration output
 *                         (\ref YAC_YAML_EMITTER_DEFAULT or
 *                          \ref YAC_YAML_EMITTER_JSON)
 *  @param[out] config     coupling configuration
 *  @remark the user is responsible for freeing the returned coupling
 *          configuration
 */
void yac_cenddef_and_emit_config_instance (
  int yac_instance_id, int emit_flags, char ** config);

/* --------------------------------------------------------------------------------
           query routines
   -------------------------------------------------------------------------------- */

/** query routine for the start datetime of the job of the default YAC instance
*/

char * yac_cget_start_datetime ( void );

/** query routine for the start datetime of the job

     @param[in] yac_instance_id id of the YAC instance
*/

char * yac_cget_start_datetime_instance ( int yac_instance_id );

/** query routine for the end datetime of the job of the default YAC instance
*/

char * yac_cget_end_datetime ( void );

/** query routine for the end datetime of the job

     @param[in] yac_instance_id id of the YAC instance
*/

char * yac_cget_end_datetime_instance ( int yac_instance_id );

/** query routine for the yac version
*/

char * yac_cget_version ( void );

/* -------------------------------------------------------------------------------- */

/** query routine for number of components defined in the default YAC
 *  instance

     @return   number of components

*/
int yac_cget_nbr_comps ( void );

/** query routine for number of components

    @param[in] yac_instance_id   id of the YAC instance
    @return                      number of components
*/
int yac_cget_nbr_comps_instance ( int yac_instance_id );

/** query routine for number of grids defined in the default YAC
 *  instance that are referenced by at least one defined field or have
 *  metadata assigned to it

     @return                  number of fields defined on the calling process

*/
int yac_cget_nbr_grids ( );

/** query routine for number of grids that are referenced by at least
 *  one defined field or have metadata assigned to it

    @param[in] yac_instance_id   id of the YAC instance
    @return                      number of fields defined on the calling process
*/
int yac_cget_nbr_grids_instance ( int yac_instance_id );

/** query routine for number of grids defined in the default YAC
 *  instance for a given components that is referenced by at least
 *  one defined field.

    @param[in] comp_name
    @return                      number of components
*/
int yac_cget_comp_nbr_grids ( const char* comp_name );

/** query routine for number of grids defined for a given components
 *  that is referenced by at least one defined field.

    @param[in] yac_instance_id   id of the YAC instance
    @param[in] comp_name
    @return                      number of components
*/
int yac_cget_comp_nbr_grids_instance ( int yac_instance_id,
  const char* comp_name );

/** query routine for number of coupling fields defined in the given component
    and grids of the default YAC instance

     @param[in] comp_name     name of the component
     @param[in] grid_name     name of the grid
     @return                  number of fields defined on the calling process

*/
int yac_cget_nbr_fields ( const char* comp_name, const char* grid_name );

/** query routine for number of coupling fields defined in the given component
    and grids

    @param[in] yac_instance_id   id of the YAC instance
    @param[in] comp_name         name of the component
    @param[in] grid_name         name of the grid
    @return                      number of fields defined on the calling process

*/
int yac_cget_nbr_fields_instance ( int yac_instance_id, const char* comp_name,
  const char* grid_name);

/* -------------------------------------------------------------------------------- */

/** query routine to get the list of component names for defined the
 *  default YAC instance.

     @param[in]  nbr_comps     number of components
     @param[out] comp_names    list of component names
*/
void yac_cget_comp_names ( int nbr_comps, const char ** comp_names );

/** query routine to get the list of component names for defined

     @param[in] yac_instance_id id of the YAC instance
     @param[in]  nbr_comps     number of components
     @param[out] comp_names  list of component names

*/
void yac_cget_comp_names_instance ( int yac_instance_id, int nbr_comps,
                                    const char ** comp_names );

/** query routine to get the list of grid names of
 *  the default YAC instance that are referenced by at least
 *  one defined field or have metadata assigned to it

     @param[in]  nbr_grids     number of grids
     @param[out] grid_names    list of grid names

*/
void yac_cget_grid_names ( int nbr_grids, const char ** grid_names );

/** query routine to get the list of grid names that are referenced by at least
 *  one defined field or have metadata assigned to it

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  nbr_grids     number of grids
     @param[out] grid_names      list of grid names

*/
void yac_cget_grid_names_instance ( int yac_instance_id,
                                    int nbr_grids, const char ** grid_names );

/** query routine to get the list of grid names of a given component
 *  that are referenced by at least one defined field or have metadata
 *  assigned to it

    @param[in]  comp_name       name of the component
    @param[in]  nbr_grids       number of grids
    @param[out] grid_names      list of grid names

*/
void yac_cget_comp_grid_names ( const char* comp_name,
                                int nbr_grids, const char ** grid_names );

/** query routine to get the list of grid names of a given component
 *  that are referenced by at least one defined field or have metadata
 *  assigned to it

    @param[in]  yac_instance_id id of the YAC instance
    @param[in]  comp_name       name of the component
    @param[in]  nbr_grids       number of grids
    @param[out] grid_names      list of grid names

*/
void yac_cget_comp_grid_names_instance ( int yac_instance_id, const char* comp_name,
                                         int nbr_grids, const char ** grid_names );

/** query routine to get the list of field names defined on the given component
    and grid for the default YAC instance.

     @param[in]  comp_name     component name
     @param[in]  grid_name     grid name
     @param[in]  nbr_fields    number of fields
     @param[out] field_names   list of field names

*/
void yac_cget_field_names ( const char* comp_name, const char* grid_name,
                            int nbr_fields, const char ** field_names );

/** query routine to get the list of field names on the given component and grid

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name       component name
     @param[in]  grid_name       grid name
     @param[in]  nbr_fields      number of fields
     @param[out] field_names      list of grid names

*/
void yac_cget_field_names_instance ( int yac_instance_id,
  const char * comp_name, const char* grid_name, int nbr_fields,
  const char ** field_names );

/* ---------------------------------------------------------------------- */

/** query routine to get the component name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return comp_name     component name of the fields component.
                               (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_component_name_from_field_id ( int field_id );

/** query routine to get the grid name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return grid_name     grid name of the fields grid.
                               (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_grid_name_from_field_id ( int field_id );

/** query routine to get the field name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return field_name    field name that has been given to yac_cdef_field.
                               (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_field_name_from_field_id ( int field_id );

/** query routine to get the timestep of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return timestep  timestep that has been given to yac_cdef_field
                           iso8601 string.
                           (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_timestep_from_field_id ( int field_id );

/** query routine to get the collection_size of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return collection_size

*/

int yac_cget_collection_size_from_field_id ( int field_id );

/** query routine to get the role of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return role          exchange type of the field\n
                           (\ref YAC_EXCHANGE_TYPE_NONE,
                            \ref YAC_EXCHANGE_TYPE_SOURCE, or
                            \ref YAC_EXCHANGE_TYPE_TARGET)
     @remark the role is set by the call to \ref yac_cenddef or
             \ref yac_cenddef_instance

*/

int yac_cget_role_from_field_id ( int field_id );

/** query routine to get the field_id from component, grid and field
    name (if defined on this process)

     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return field_id
*/
int yac_cget_field_id(const char* comp_name, const char* grid_name,
                      const char* field_name);

/** query routine to get the field_id from component, grid and field
    name (if defined on this process)

     @param[in] yac_instance_id
     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return field_id
*/
int yac_cget_field_id_instance(int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name);

/* ---------------------------------------------------------------------- */

/** Get metadata for a component
    @param[in] comp_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_component_metadata(const char* comp_name);

/** Get metadata for a component
    @param[in] yac_instance_id
    @param[in] comp_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_component_metadata_instance(int yac_instance_id,
  const char* comp_name);

/** Get metadata for a grid
    @param[in] grid_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_grid_metadata(const char* grid_name);

/** Get metadata for a grid
    @param[in] yac_instance_id
    @param[in] grid_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_grid_metadata_instance(int yac_instance_id,
  const char* grid_name);

/** Get metadata for a field
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_field_metadata(const char* comp_name, const char* grid_name,
  const char* field_name);

/** Get metadata for a field
    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return metadata (NULL if no metadata is defined)
*/
const char* yac_cget_field_metadata_instance(int yac_instance_id,
  const char* comp_name, const char* grid_name, const char* field_name);

/* ---------------------------------------------------------------------- */

/** query routine to get the timestep of a field

    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return timestep timestep that has been given to yac_cdef_field iso8601 string
                         (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_field_timestep ( const char* comp_name, const char* grid_name,
  const char* field_name );

/** query routine to get the timestep of a field

    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return timestep timestep that has been given to yac_cdef_field iso8601 string
                         (array of length YAC_MAX_CHARLEN)

*/

const char* yac_cget_field_timestep_instance ( int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name );

/** query routine to get the fractional mask fallback value of a field

    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return frac_mask_fallback_value
    @remark if YAC_FRAC_MASK_NO_VALUE is returned, then no
            fractional mask fallback value is defined for this field

*/

double yac_cget_field_frac_mask_fallback_value(
  const char* comp_name, const char* grid_name, const char* field_name);

/** query routine to get the collection_size of a field

    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return collection_size

*/

int yac_cget_field_collection_size(
  const char* comp_name, const char* grid_name, const char* field_name);

/** query routine to get the fractional mask fallback value of a field

    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return frac_mask_fallback_value
    @remark if \ref YAC_FRAC_MASK_NO_VALUE is returned, then no
            fractional mask fallback value is defined for this field

*/

double yac_cget_field_frac_mask_fallback_value_instance(
  int yac_instance_id, const char* comp_name, const char* grid_name,
  const char* field_name );

/** query routine to get the collection_size of a field

    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return collection_size

*/

int yac_cget_field_collection_size_instance(
  int yac_instance_id, const char* comp_name, const char* grid_name,
  const char* field_name);

/** query routine to get the role of a field for a given ID

    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return role          exchange type of the field\n
                          (\ref YAC_EXCHANGE_TYPE_NONE,
                          \ref YAC_EXCHANGE_TYPE_SOURCE, or
                          \ref YAC_EXCHANGE_TYPE_TARGET)
    @remark the role is set by the call to \ref yac_cenddef or
                          \ref yac_cenddef_instance

*/

int yac_cget_field_role ( const char* comp_name, const char* grid_name,
  const char* field_name );

/** query routine to get the role of a field for a given ID

    @param[in] yac_instance_id
    @param[in] comp_name
    @param[in] grid_name
    @param[in] field_name
    @return role          exchange type of the field\n
                          (\ref YAC_EXCHANGE_TYPE_NONE,
                          \ref YAC_EXCHANGE_TYPE_SOURCE, or
                          \ref YAC_EXCHANGE_TYPE_TARGET)
    @remark the role is set by the call to \ref yac_cenddef or
                          \ref yac_cenddef_instance

*/

int yac_cget_field_role_instance ( int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name );

/* --------------------------------------------------------------------------------
           auxiliary routines
   -------------------------------------------------------------------------------- */

#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif

/**
 * functions used as error handler must conform to this interface
 */
typedef void (*yac_abort_func)(MPI_Comm comm, const char *msg,
                               const char *source, int line)
  __attribute__((noreturn));

/**
 * redirect stderr and stdout into files with names depending on the process ID
 *
 * The original version was part of the OASIS4 coupler and has been adjusted to
 * work for yac.
 *
 * @param[in]  filestem   basename of stdout/stderr files
 * @param[in]  parallel   flag for direction of output
                             parallel = 1 : Create single output files for each process
                             otherwise    : Create output files only for process 0.
                                           The output of other processes is not redirected.
 * @param[in]  my_pe      local rank
 * @param[in]  npes       total number of processes in set
 */

void yac_redirstdout ( const char *filestem,
                       int parallel, int my_pe,
                       int npes );

/**
 * Calls the currently set abort handler (\ref abort_default "yac_abort_default" by default)
 * @param[in] comm   MPI communicator used to call MPI_Abort
 * @param[in] msg    message text to print
 * @param[in] source string describing source file name
 * @param[in] line   line number of caller
 */
void yac_abort(MPI_Comm comm, const char *msg,
               const char *source, int line)
  __attribute__((noreturn));

/**
 * Call the \ref yac_abort function (providing the default communicator for
 * the comm argument).
 * @param msg    message text to print
 * @param source string describing source file name
 * @param line   line number of caller
 */
void yac_abort_message(char const *msg, const char *source, int line);

/**
 * Restores default abort handler
 */
void yac_restore_default_abort_handler(void);

/**
 * Sets custom abort handler
 * @param[in] custom_abort custom abort handler
 * @remark This abort handler should call MPI_Abort if possible. Once this
 *         handler was called, the internal state of YAC is undefined and no
 *         YAC functions should be called anymore.
 */
void yac_set_abort_handler(yac_abort_func custom_abort);

/**
 * Gets abort handler
 * @return currently set abort handler
 */
yac_abort_func yac_get_abort_handler(void);

/**
 * Gets default abort handler
 * @return default abort handler
 */
yac_abort_func yac_get_default_abort_handler(void);

/**
 * Sets default MPI communicator (MPI_COMM_WORLD by default)
 * @param[in] comm default MPI communicator
 */
void yac_set_default_comm(MPI_Comm comm);

/* -------------------------------------------------------------------------------- */

/** query routine to get the number of points in a grid

     @param[in] location
     @param[in] grid_id

*/

size_t yac_get_grid_size ( int location, int grid_id );

/* -------------------------------------------------------------------------------- */

/** query routine to get the number of points in a pointset

     @param[in] points_id

*/

size_t yac_get_points_size ( int points_id  );

/* -------------------------------------------------------------------------------- */

/** gets an empty stack trace
     @param[out] interp_stack_config_id interpolation stack
*/
void yac_cget_interp_stack_config(int * interp_stack_config_id);

/** frees a stack trace
     @param[in] interp_stack_config_id interpolation stack
*/
void yac_cfree_interp_stack_config(int interp_stack_config_id);

/** adds average interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] reduction_type         reduction type
     @param[in] partial_coverage       allow partial coverage
*/
void yac_cadd_interp_stack_config_average(
  int interp_stack_config_id, int reduction_type, int partial_coverage);

/** adds n-nearest-neighbour interpolation to the buttom of a
    interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] type                   reduction type
     @param[in] n                      number of nearest neighbour points
     @param[in] scale                  scale parameter required by some
                                       reduction types
*/
void yac_cadd_interp_stack_config_nnn(
  int interp_stack_config_id, int type, size_t n, double scale);

/** adds conservative interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] order                  first or second order
     @param[in] enforced_conserv       enforce local conservation
     @param[in] partial_coverage       allow partial coverage
     @param[in] normalisation          normalisation type
*/
void yac_cadd_interp_stack_config_conservative(
  int interp_stack_config_id, int order, int enforced_conserv,
  int partial_coverage, int normalisation);

/** adds source point mapping interpolation to the buttom
    of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] spread_distance        spread distance
     @param[in] max_search_distance    maximum search distance
     @param[in] weight_type            reduction type
*/
void yac_cadd_interp_stack_config_spmap(
  int interp_stack_config_id, double spread_distance,
  double max_search_distance, int weight_type);

/** adds HCSBB interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
*/
void yac_cadd_interp_stack_config_hcsbb(int interp_stack_config_id);

/** adds user weight file interpolation to the buttom
    of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] filename               weight file name
     @param[in] src_grid_name          name of source grid
     @param[in] tgt_grid_name          name of target grid
     @remark YAC will read the wight file in parallel.
             This parallel output can be configured as described here:
             \ref io_config_detail
*/
void yac_cadd_interp_stack_config_user_file(
  int interp_stack_config_id, char const * filename, char const * src_grid_name,
  char const * tgt_grid_name);

/** adds fixed interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] value                  fixed value
*/
void yac_cadd_interp_stack_config_fixed(
  int interp_stack_config_id, double value);

/** adds fixed interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] constructor_key        key provided to
                                       \ref yac_interp_method_check_add_constructor_callback
                                       for a constructor callback routine
     @param[in] do_search_key          key provided to
                                       \ref yac_interp_method_check_add_do_search_callback
                                       for a do_search callback routine
*/
void yac_cadd_interp_stack_config_check(
  int interp_stack_config_id, char const * constructor_key,
  char const * do_search_key);

/** adds creep interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] creep_distance         creep distance
*/
void yac_cadd_interp_stack_config_creep(
  int interp_stack_config_id, int creep_distance);

/** adds user callback interpolation to the buttom of a interpolation stack
     @param[in] interp_stack_config_id   interpolation stack
     @param[in] func_compute_weights_key key provided to
                                         \ref yac_cadd_compute_weights_callback
                                         for a compute_weights callback routine
*/
void yac_cadd_interp_stack_config_user_callback(
  int interp_stack_config_id, char const * func_compute_weights_key);

/** method signature for weight compuation function used by
    \ref interp_method_callback

  If the field value of a target point, whose weights have been computed
  by this routine, is to be computed on the source process at which this
  routine was called, YAC would use the weights as follows (it is assumed
  that all required source points are available at this process):

  \code{.c}

  // this routine return an index to a source field value associated to
  // the provided global id
  size_t global_id_2_local_idx(int src_field_idx, int global_id);

  [...]

  double tgt_value[collection_size];
  for (int collection_idx = 0; collection_idx < collection_size;
       collection_idx)
    tgt_value[collection_idx] = 0.0;

  // for all source fields
  for (int src_field_idx = 0; src_field_idx < num_pointsets; ++src_field_idx) {

    for (size_t i = 0; i < result_count[src_field_idx]; ++i) {

      size_t src_field_point_idx =
        global_id_2_local_idx(
          src_field_idx, global_results_points[src_field_idx][i]);
      double weight = result_weights[src_field_idx][src_field_point_idx];

      for (int collection_idx = 0; collection_idx < collection_size;
           collection_idx) {

        tgt_value[collection_idx] +=
          src_field[collection_idx][src_field_idx][src_field_point_idx] *
          weight;
      }
    }
  }

  \endcode

  If one or more source points associated to the global ids return by this
  routine are not actually available at this process, YAC will ensure that
  the respective value will be available for the interpolation.

  @param[in]  tgt_coords            3D coordinates of the target point
  @param[in]  src_cell_id           global id of the source cell matching
                                    the target point
  @param[in]  src_cell_idx          (zero based) index of the source
                                    cell matching the target point
  @param[out] global_results_points global ids of source points to be used
                                    for the interpolation of the target point
  @param[out] result_weights        weights to be used for the interpolation
                                    of the target point
  @param[out] result_count          number of source points to be used
                                    for the interpolation of the target point
                                    per pointset
  @param[in]  user_data             user_data pointer provided to
                                    \ref yac_cadd_compute_weights_callback
  @remark The size of the first dimension of the arrays global_results_points,
          result_weights, and result_count has to be equal to num_pointsets
          provided to \ref yac_cdef_field or \ref yac_cdef_field_mask that was
          used to define the field involved in the coupling that is
          supposed to use this weight computation function.
  @remark 3D coordinates are for a point of a unit sphere
*/
typedef void (*yac_func_compute_weights)(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);

/** registers a callback routine for the computation of weights by the
 *  interpolation method user_callback

     @param[in] compute_weights_callback pointer to a weight computation routine
     @param[in] user_data                data pointer that will be passed to
                                         compute_weights_callback
     @param[in] key                      key for identifing the callback routine
     @remark the callback has to be set on all source processes taking part in
             respective coupling
     @remark the key has to match the one provided in coupling configuration
*/
void yac_cadd_compute_weights_callback(
  yac_func_compute_weights compute_weights_callback,
  void * user_data, char const * key);

#endif // YAC_INTERFACE_H

