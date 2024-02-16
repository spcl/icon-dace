/**
 * @file yac_interface.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Teresa Holfeld <teresa.holfeld@zmaw.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Teresa Holfeld <teresa.holfeld@zmaw.de>
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ensure_array_size.h"
#include "event.h"
#include "interpolation.h"
#include "interp_method_callback.h"
#include "yac_interface.h"
#include "utils.h"
#include "version.h"
#include "interp_grid.h"
#include "instance.h"
#include "fields.h"
#include "couple_config.h"
#include "config_yaml.h"

int const YAC_LOCATION_CELL   = CELL;
int const YAC_LOCATION_CORNER = CORNER;
int const YAC_LOCATION_EDGE   = EDGE;

int const YAC_EXCHANGE_TYPE_NONE   = NOTHING;
int const YAC_EXCHANGE_TYPE_SOURCE = SOURCE;
int const YAC_EXCHANGE_TYPE_TARGET = TARGET;

int const YAC_ACTION_NONE               = NONE;
int const YAC_ACTION_REDUCTION          = REDUCTION;
int const YAC_ACTION_COUPLING           = COUPLING;
int const YAC_ACTION_RESTART            = RESTART;
int const YAC_ACTION_GET_FOR_RESTART    = GET_FOR_RESTART;
int const YAC_ACTION_PUT_FOR_RESTART    = PUT_FOR_RESTART;
int const YAC_ACTION_GET_FOR_CHECKPOINT = GET_FOR_CHECKPOINT;
int const YAC_ACTION_PUT_FOR_CHECKPOINT = PUT_FOR_CHECKPOINT;
int const YAC_ACTION_OUT_OF_BOUND       = OUT_OF_BOUND;

int const YAC_REDUCTION_TIME_NONE       = TIME_NONE;
int const YAC_REDUCTION_TIME_ACCUMULATE = TIME_ACCUMULATE;
int const YAC_REDUCTION_TIME_AVERAGE    = TIME_AVERAGE;
int const YAC_REDUCTION_TIME_MINIMUM    = TIME_MINIMUM;
int const YAC_REDUCTION_TIME_MAXIMUM    = TIME_MAXIMUM;

int const YAC_TIME_UNIT_MILLISECOND = C_MILLISECOND;
int const YAC_TIME_UNIT_SECOND = C_SECOND;
int const YAC_TIME_UNIT_MINUTE = C_MINUTE;
int const YAC_TIME_UNIT_HOUR = C_HOUR;
int const YAC_TIME_UNIT_DAY = C_DAY;
int const YAC_TIME_UNIT_MONTH = C_MONTH;
int const YAC_TIME_UNIT_YEAR = C_YEAR;
int const YAC_TIME_UNIT_ISO_FORMAT = C_ISO_FORMAT;

int const YAC_CALENDAR_NOT_SET = CALENDAR_NOT_SET;
int const YAC_PROLEPTIC_GREGORIAN = PROLEPTIC_GREGORIAN;
int const YAC_YEAR_OF_365_DAYS = YEAR_OF_365_DAYS;
int const YAC_YEAR_OF_360_DAYS = YEAR_OF_360_DAYS;

int const YAC_AVG_ARITHMETIC = AVG_ARITHMETIC;
int const YAC_AVG_DIST = AVG_DIST;
int const YAC_AVG_BARY = AVG_BARY;

int const YAC_NNN_AVG = NNN_AVG;
int const YAC_NNN_DIST = NNN_DIST;
int const YAC_NNN_GAUSS = NNN_GAUSS;
int const YAC_NNN_RBF = NNN_RBF;

int const YAC_CONSERV_DESTAREA = CONSERV_DESTAREA;
int const YAC_CONSERV_FRACAREA = CONSERV_FRACAREA;

int const YAC_SPMAP_AVG = SPMAP_AVG;
int const YAC_SPMAP_DIST = SPMAP_DIST;

struct user_input_data_component {
  struct yac_instance * instance;
  size_t component_index;
  char * name;
};

struct user_input_data_points {
  int default_mask_id;
  struct yac_basic_grid * grid;
  enum yac_location location;
  size_t coordinates_idx;
};

struct user_input_data_masks {
  struct yac_basic_grid * grid;
  enum yac_location location;
  size_t masks_idx;
};

static int default_instance_id = INT_MAX;

static int yac_instance_count = 0;

void ** pointer_lookup_table = NULL;
static int pointer_lookup_table_size = 0;

struct yac_basic_grid ** grids = NULL;
size_t num_grids = 0;

static struct user_input_data_component ** components = NULL;
static size_t num_components = 0;

static struct user_input_data_points ** points = NULL;
static size_t num_points = 0;

static struct user_input_data_masks ** masks = NULL;
static size_t num_masks = 0;

static struct yac_interp_stack_config ** interp_stack_configs = NULL;
static size_t num_interp_stack_configs = 0;

struct yac_ext_couple_config {
  char * weight_file;
  int mapping_side;
  double scale_factor;
  double scale_summand;
  size_t num_src_mask_names;
  char ** src_mask_names;
  char * tgt_mask_name;
};

/**
 * gives a unique index for a given pointer
 * @param[in] pointer
 * @return unique value associated to pointer
 * @see unique_id_to_pointer
 */
static int yac_pointer_to_unique_id(void * pointer) {

  pointer_lookup_table =
    xrealloc (pointer_lookup_table,
      (pointer_lookup_table_size + 1) * sizeof(*pointer_lookup_table));

  pointer_lookup_table[pointer_lookup_table_size] = pointer;

  pointer_lookup_table_size++;

  return pointer_lookup_table_size - 1;
}

/**
 * gives the pointer that is associated to the given id
 * @param[in] id      unique index previously returned by pointer_to_unique_id
 * @param[in] id_name name of the id
 * @return pointer the is associated to the given id \n NULL if the id is invalid
 */
static void * yac_unique_id_to_pointer(int id, char const * id_name) {

  YAC_ASSERT_F(
    id < pointer_lookup_table_size &&
    id >= 0,
    "ERROR(yac_unique_id_to_pointer): invalid %s", id_name)
  return pointer_lookup_table[id];
}

/**
 * gives the unique index associated with pointer
 * @param[in] pointer
 * @return unique value associated to pointer
 * @remark the pointer has to have been registers using
 *         \ref yac_pointer_to_unique_id, otherwise this routine returns
 *         INT_MAX
 */
static int yac_lookup_pointer(void const * pointer) {

  int ret = INT_MAX;
  for (int i = 0; (ret == INT_MAX) && (i < pointer_lookup_table_size); ++i)
    if (pointer_lookup_table[i] == pointer) ret = i;
  return ret;
}

/**
 * frees all memory used for the pointer/unique_id conversion
 * \remarks this should only be called after the last call to
 *          \ref yac_pointer_to_unique_id and \ref yac_unique_id_to_pointer, because afterwards
 *          \ref yac_unique_id_to_pointer will not be able to return the respective pointers
 *          for previously valid unique ids
 */
static void yac_free_pointer_unique_lookup() {

  free(pointer_lookup_table);
  pointer_lookup_table = NULL;
  pointer_lookup_table_size = 0;
}

/* ---------------------------------------------------------------------- */

static int yac_add_grid(
  const char * grid_name, struct basic_grid_data grid_data) {

  for (size_t i = 0; i < num_grids; ++i)
    YAC_ASSERT(
      strcmp(yac_basic_grid_get_name(grids[i]), grid_name),
      "ERROR(yac_add_grid): multiple definitions of grid with identical name")

  YAC_ASSERT(
    strlen(grid_name) <= YAC_MAX_CHARLEN,
    "ERROR(yac_add_grid): grid name is too long (maximum is YAC_MAX_CHARLEN)")

  struct yac_basic_grid * grid = yac_basic_grid_new(grid_name, grid_data);

  grids = xrealloc(grids, (num_grids + 1) * sizeof(*grids));
  grids[num_grids] = grid;
  num_grids++;
  return yac_pointer_to_unique_id(grid);
}

static void yac_free_grids() {

  for (size_t i = 0; i < num_grids; ++i)
    yac_basic_grid_delete(grids[i]);
  free(grids);
  grids = NULL;
  num_grids = 0;
}

/* ---------------------------------------------------------------------- */

static void yac_free_points() {

  for (size_t i = 0; i < num_points; ++i) free(points[i]);
  free(points);
  points = NULL;
  num_points = 0;
}

/* ---------------------------------------------------------------------- */

static void yac_free_masks() {

  for (size_t i = 0; i < num_masks; ++i) free(masks[i]);
  free(masks);
  masks = NULL;
  num_masks = 0;
}

/* ---------------------------------------------------------------------- */

static void yac_free_components() {

  for (size_t i = 0; i < num_components; ++i){
    free(components[i]->name);
    free(components[i]);
  }
  free(components);
  components = NULL;
  num_components = 0;
}

/* ---------------------------------------------------------------------- */

static void yac_free_interp_stack_configs() {

  for (size_t i = 0; i < num_interp_stack_configs; ++i)
    if (interp_stack_configs[i] != NULL)
      yac_interp_stack_config_delete(interp_stack_configs[i]);
  free(interp_stack_configs);
  interp_stack_configs = NULL;
  num_interp_stack_configs = 0;
}

/* ---------------------------------------------------------------------- */

static void check_default_instance_id(char const * routine_name) {

  YAC_ASSERT_F(
    default_instance_id != INT_MAX,
    "ERROR(%s): no default yac instance is defined yet", routine_name)
}

/* ---------------------------------------------------------------------- */

static void yac_check_version(MPI_Comm comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  char recv_version[64];

  size_t yac_version_len = strlen(yac_version) + 1;

  YAC_ASSERT_F(
    yac_version_len <= sizeof(recv_version),
    "ERROR(yac_check_version): "
    "version string \"%s\" is too long (has to be shorter than %zu)",
    yac_version, sizeof(recv_version))

  if (comm_rank == 0) strcpy(recv_version, yac_version);

  yac_mpi_call(
    MPI_Bcast(recv_version, (int)yac_version_len, MPI_CHAR, 0, comm), comm);

  YAC_ASSERT_F(
    !strcmp(recv_version, yac_version),
    "ERROR(yac_check_version): inconsistent yac versions between processes "
    "(on local process \"%s\"; on root \"%s\")", yac_version, recv_version)
}

void check_mpi_initialised(char const * routine_name) {

  YAC_ASSERT_F(
    yac_mpi_is_initialised(),
    "ERROR(%s): "
    "MPI has not yet been initialised, but this call got a "
    "communicator which is not allowed by the MPI standard "
    "(MPI_COMM_WORLD is not to be used before the call to "
    "MPI_Init or MPI_Init_thread)", routine_name);
}

void yac_cmpi_handshake(MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms) {
  check_mpi_initialised("yac_cmpi_handshake");
  yac_mpi_handshake(comm, n, group_names, group_comms);
}

void yac_cmpi_handshake_f2c(MPI_Fint comm, int n, char const** group_names,
  MPI_Fint * group_comms) {
  MPI_Comm comm_c = MPI_Comm_f2c(comm);
  MPI_Comm * group_comms_c = xmalloc(n*sizeof(*group_comms_c));
  yac_mpi_handshake(comm_c, n, group_names, group_comms_c);
  for(int i = 0; i<n; ++i)
    group_comms[i] = MPI_Comm_c2f(group_comms_c[i]);
  free(group_comms_c);
}

static int yac_init(MPI_Comm comm) {

  check_mpi_initialised("yac_init");
  yac_yaxt_init(comm);
  yac_check_version(comm);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  yac_instance_count++;
  return
    yac_pointer_to_unique_id(
      yac_instance_new(comm));
}

int yac_cget_default_instance_id(){
  check_default_instance_id("yac_cget_default_instance_id");
  return default_instance_id;
}

void yac_cinit_comm_instance (
  MPI_Comm comm, int * yac_instance_id) {

  *yac_instance_id = yac_init(comm);
}

void yac_cinit_comm_instance_f2c (
  MPI_Fint comm, int * yac_instance_id) {

  yac_cinit_comm_instance(MPI_Comm_f2c(comm), yac_instance_id);
}

void yac_cinit_comm ( MPI_Comm comm ) {

  YAC_ASSERT(
    default_instance_id == INT_MAX,
    "ERROR(yac_cinit_comm): default yac instance already defined")

  yac_cinit_comm_instance(comm, &default_instance_id);
}

void yac_cinit_comm_f2c ( MPI_Fint comm_f ) {

  yac_cinit_comm(MPI_Comm_f2c(comm_f));
}

void yac_cinit_instance(
  int * yac_instance_id) {

  // we have to initalise MPI here in order to be able to use MPI_COMM_WORLD
  yac_mpi_init();
  MPI_Comm yac_comm;
  const char* group_name = "yac";
  yac_cmpi_handshake(MPI_COMM_WORLD, 1, &group_name, &yac_comm);
  yac_cinit_comm_instance(yac_comm, yac_instance_id);
  yac_mpi_call(MPI_Comm_free(&yac_comm), MPI_COMM_WORLD);
}

void yac_cinit ( void ) {
  YAC_ASSERT(
    default_instance_id == INT_MAX,
    "ERROR(yac_cinit_comm): default yac instance already defined")

  yac_cinit_instance(&default_instance_id);
}

/* ---------------------------------------------------------------------- */

void yac_cinit_comm_dummy ( MPI_Comm comm ) {

  check_mpi_initialised("yac_cinit_comm_dummy");

  yac_yaxt_init(comm);
  yac_check_version(comm);

  yac_instance_dummy_new(comm);
}

void yac_cinit_comm_dummy_f2c ( MPI_Fint comm_f ) {

  check_mpi_initialised("yac_cinit_comm_dummy_f2c");
  yac_cinit_comm_dummy(MPI_Comm_f2c(comm_f));
}

void yac_cinit_dummy( void ) {

  // we have to initalise MPI here in order to be able to use MPI_COMM_WORLD
  yac_mpi_init();
  yac_cmpi_handshake(MPI_COMM_WORLD, 0, NULL, NULL);
}

/* ---------------------------------------------------------------------- */

int yac_cyaml_get_emitter_flag_default_c2f() {
  return YAC_YAML_EMITTER_DEFAULT;
}

int yac_cyaml_get_emitter_flag_json_c2f() {
  return YAC_YAML_EMITTER_JSON;
}

/* ---------------------------------------------------------------------- */


void yac_cread_config_yaml_instance( int yac_instance_id,
  const char * yaml_filename){
  struct yac_instance * instance
    = yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  yac_yaml_read_coupling(
    yac_instance_get_couple_config(instance), yaml_filename,
    YAC_YAML_PARSER_DEFAULT);
}

void yac_cread_config_yaml( const char * yaml_filename){
  check_default_instance_id("read_config_yaml");
  yac_cread_config_yaml_instance(default_instance_id, yaml_filename);
}

void yac_cread_config_json_instance( int yac_instance_id,
  const char * yaml_filename){
  struct yac_instance * instance
    = yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  yac_yaml_read_coupling(
    yac_instance_get_couple_config(instance), yaml_filename,
    YAC_YAML_PARSER_JSON_FORCE);
}

void yac_cread_config_json( const char * yaml_filename){
  check_default_instance_id("read_config_yaml");
  yac_cread_config_json_instance(default_instance_id, yaml_filename);
}

/* ---------------------------------------------------------------------- */

static void yac_ccleanup_instance_(int yac_instance_id) {

  struct yac_instance * instance
    = yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");

  /* free memory */
  yac_instance_delete(instance);
  --yac_instance_count;
  if(yac_instance_id == default_instance_id)
    default_instance_id = INT_MAX;

  // remove dangeling pointers from components
  for(int i = 0; i<num_components;++i){
    if(components[i]->instance == instance){
      components[i]->instance = NULL;
    }
  }
}

void yac_ccleanup_instance (int yac_instance_id) {

  yac_ccleanup_instance_(yac_instance_id);

  /* cleanup mpi */
  yac_mpi_cleanup();
}

void yac_ccleanup () {

  check_default_instance_id("yac_ccleanup");

  yac_ccleanup_instance(default_instance_id);
}

/* ---------------------------------------------------------------------- */

static void cleanup() {
  yac_mpi_finalize();
  yac_free_pointer_unique_lookup();
  yac_free_grids();
  yac_free_points();
  yac_free_masks();
  yac_free_components();
  yac_free_interp_stack_configs();
}

void yac_cfinalize_instance(int yac_instance_id) {

  yac_ccleanup_instance_(yac_instance_id);

  /* cleanup close MPI */

  if (yac_instance_count == 0) cleanup();
}

void yac_cfinalize () {

  if (default_instance_id != INT_MAX)
    yac_ccleanup_instance_(default_instance_id);

  /* cleanup close MPI */

  if (yac_instance_count == 0) cleanup();
}

/* ---------------------------------------------------------------------- */

void yac_cdef_datetime_instance (
  int yac_instance_id, const char * start_datetime,
  const char * end_datetime) {

  yac_instance_def_datetime(
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"),
    start_datetime, end_datetime);
}

void yac_cdef_datetime ( const char * start_datetime,
                         const char * end_datetime ) {

  check_default_instance_id("yac_cdef_datetime");
  yac_cdef_datetime_instance(
    default_instance_id, start_datetime, end_datetime);
}

void yac_cdef_calendar ( int calendar ) {

  YAC_ASSERT(
    (calendar == YAC_PROLEPTIC_GREGORIAN) ||
    (calendar == YAC_YEAR_OF_360_DAYS) ||
    (calendar == YAC_YEAR_OF_365_DAYS),
    "ERROR(yac_cdef_calendar): invalid calendar type")

  calendarType curr_calendar = getCalendarType();
  YAC_ASSERT(
    (curr_calendar == CALENDAR_NOT_SET) ||
    (curr_calendar == (calendarType)calendar),
    "ERROR(yac_cdef_calendar): inconsistent calendar definition")

  initCalendar((calendarType)calendar);
}

/* ---------------------------------------------------------------------- */

char * yac_cget_start_datetime_instance(int yac_instance_id) {

  return
    yac_instance_get_start_datetime(
      yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"));
}

char * yac_cget_start_datetime ( void ) {

  check_default_instance_id("yac_cget_start_datetime");
  return yac_cget_start_datetime_instance(default_instance_id);
}

/* ---------------------------------------------------------------------- */

char * yac_cget_end_datetime_instance(int yac_instance_id) {

  return
    yac_instance_get_end_datetime(
      yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"));
}

char * yac_cget_end_datetime ( void ) {

  check_default_instance_id("yac_cget_end_datetime");
  return yac_cget_end_datetime_instance(default_instance_id);
}

/* ---------------------------------------------------------------------- */

/* The YAC version number is provided in version.h
 */

char * yac_cget_version ( void ) {
  return yac_version;
}

/* ---------------------------------------------------------------------- */

static struct user_input_data_component * get_user_input_data_component(
  int comp_id, char const * routine) {

  struct user_input_data_component * comp_info =
    yac_unique_id_to_pointer(comp_id, "comp_id");

  YAC_ASSERT_F(
    comp_info->instance != NULL,
    "ERROR(%s): instance of component \"%s\" is already finalized",
    routine, comp_info->name);

  YAC_ASSERT_F(
    comp_info->component_index != SIZE_MAX,
    "ERROR(%s): component definition has not yet been completed "
    "(component name: \"%s\")",
    routine, comp_info->name);

  return comp_info;
}

/* ---------------------------------------------------------------------- */

/* c interface routine */

void yac_cget_comp_comm ( int comp_id, MPI_Comm *comp_comm ) {

  struct user_input_data_component * comp_info =
    get_user_input_data_component(comp_id, "yac_cget_comp_comm");

  MPI_Comm_dup(
    yac_component_get_communicator(
      yac_instance_get_component(
        comp_info->instance, comp_info->component_index)), comp_comm);
}

/* internal routine to serve the Fortran interface */

void yac_get_comp_comm_f2c ( int comp_id, MPI_Fint *comp_comm_f ) {

  MPI_Comm comp_comm;
  yac_cget_comp_comm(comp_id, &comp_comm);

  *comp_comm_f = MPI_Comm_c2f(comp_comm);
}

/* ---------------------------------------------------------------------- */

/* c interface routine */

void yac_cget_comps_comm_instance(
  int yac_instance_id,
  char const ** comp_names, int num_comps, MPI_Comm * comps_comm) {

  *comps_comm =
    yac_instance_get_comps_comm(
      yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"),
      comp_names, (size_t)num_comps);
}

void yac_cget_comps_comm(
  const char ** comp_names, int num_comps, MPI_Comm * comps_comm) {

  check_default_instance_id("yac_cget_comps_comm");
  yac_cget_comps_comm_instance(
    default_instance_id, comp_names, num_comps, comps_comm);
}

/* internal routine to serve the Fortran interface */

void yac_cget_comps_comm_instance_f2c(
  int yac_instance_id,
  char const ** comp_names, int num_comps, MPI_Fint * comps_comm_f) {

  MPI_Comm comps_comm;
  yac_cget_comps_comm_instance(
    yac_instance_id, comp_names, num_comps, &comps_comm);
  *comps_comm_f = MPI_Comm_c2f(comps_comm);
}

void yac_cget_comps_comm_f2c(
  char const ** comp_names, int num_comps, MPI_Fint *comps_comm_f) {

  MPI_Comm comps_comm;
  yac_cget_comps_comm(comp_names, num_comps, &comps_comm);
  *comps_comm_f = MPI_Comm_c2f(comps_comm);
}

/* ---------------------------------------------------------------------- */

void yac_cpredef_comp_instance(
  int yac_instance_id, char const * name, int * comp_id){

  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");

  YAC_ASSERT(
    !yac_instance_components_are_defined(instance),
    "ERROR(yac_cpredef_comp_instance): components have already been defined");

  YAC_ASSERT(
    (name != NULL) && (*name != '\0'),
    "ERROR(yac_cpredef_comp_instance): missing component name");

  for (size_t i = 0; i < num_components; ++i) {
    YAC_ASSERT_F(
      (components[i]->instance != instance) ||
      strcmp(components[i]->name, name),
      "ERROR(yac_cpredef_comp_instance): "
      "component \"%s\" is already defined", name);
  }

  components =
    xrealloc(
      components, (num_components + 1) * sizeof(*components));

  components[num_components] = xmalloc(1 * sizeof(**components));
  components[num_components]->instance = instance;
  components[num_components]->component_index = SIZE_MAX; // not determined yet
  components[num_components]->name = strdup(name);
  *comp_id = yac_pointer_to_unique_id(components[num_components]);
  num_components++;
}

void yac_cpredef_comp(char const * name, int * comp_id){
  check_default_instance_id("yac_cpredef_comp");
  yac_cpredef_comp_instance(default_instance_id, name, comp_id);
}

void yac_cdef_comps_instance(
  int yac_instance_id, char const ** comp_names, int num_comps,
  int * comp_ids) {

  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");

  for(int i = 0; i<num_comps; ++i)
    yac_cpredef_comp_instance(yac_instance_id, comp_names[i], comp_ids + i);

  // count predefined components and set their index
  size_t comp_counter = 0;
  char const ** all_comp_names = xmalloc(num_components * sizeof(*all_comp_names));
  for(int i = 0; i<num_components; ++i){
    if(components[i]->instance == instance){
      all_comp_names[comp_counter] = components[i]->name;
      components[i]->component_index = comp_counter;
      comp_counter++;
    }
  }

  yac_instance_def_components(instance, all_comp_names, comp_counter);

  free(all_comp_names);
}

void yac_cdef_comps (
  char const ** comp_names, int num_comps, int * comp_ids ) {

  check_default_instance_id("yac_cdef_comps");
  yac_cdef_comps_instance(
    default_instance_id, comp_names, num_comps, comp_ids);
}

void yac_cdef_comp_instance(
  int yac_instance_id, char const * comp_name, int * comp_id ) {

  yac_cdef_comps_instance(yac_instance_id, &comp_name, 1, comp_id);
}

void yac_cdef_comp ( char const * comp_name, int * comp_id ) {

  yac_cdef_comps(&comp_name, 1, comp_id);
}

/* ---------------------------------------------------------------------- */

static int user_input_data_add_points(
  struct yac_basic_grid * grid, enum yac_location location,
  coordinate_pointer coordinates) {

  struct user_input_data_points * curr_points = xmalloc(1 * sizeof(*curr_points));
  curr_points->default_mask_id = INT_MAX;
  curr_points->grid = grid;
  curr_points->location = location;
  curr_points->coordinates_idx =
    yac_basic_grid_add_coordinates_nocpy(grid, location, coordinates);

  points = xrealloc(points, (num_points + 1) * sizeof(*points));
  points[num_points] = curr_points;
  num_points++;

  return yac_pointer_to_unique_id(curr_points);
}

/* ---------------------------------------------------------------------- */

static void check_x_vertices(
  double const * x_vertices, size_t count, char const * routine_name) {

  for (size_t i = 0; i < count; ++i)
    YAC_ASSERT_F(
      (x_vertices[i] >= -720.0) && (x_vertices[i] <= 720.0),
      "ERROR(%s): x_vertices[%zu] = %lf outside of valid range [-720;720]",
      routine_name, i, x_vertices[i]);
}

static void check_y_vertices(
  double const * y_vertices, size_t count, char const * routine_name) {

  for (size_t i = 0; i < count; ++i)
    YAC_ASSERT_F(
      (y_vertices[i] >= -90.00001) && (y_vertices[i] <= 90.00001),
      "ERROR(%s): y_vertices[%zu] = %lf "
      "outside of valid range [90.00001;90.00001]",
      routine_name, i, y_vertices[i]);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_points_reg2d ( int const grid_id,
                             int const *nbr_points,
                             int const located,
                             double const *x_points,
                             double const *y_points,
                             int *point_id ) {

  enum yac_location location = yac_get_location(located);
  struct yac_basic_grid * grid = yac_unique_id_to_pointer(grid_id, "grid_id");

  size_t num_points = yac_basic_grid_get_data_size(grid, location);
  YAC_ASSERT(
    ((size_t)nbr_points[0] * (size_t)nbr_points[1]) == num_points,
    "ERROR(yac_cdef_points_reg2d): nbr_points does not match with grid")

  check_x_vertices(
    x_points, (size_t)(nbr_points[0]), "yac_cdef_points_reg2d");
  check_y_vertices(
    y_points, (size_t)(nbr_points[1]), "yac_cdef_points_reg2d");

  coordinate_pointer coordinates =
    xmalloc(num_points * sizeof(*coordinates));
  for (int i = 0, k = 0; i < nbr_points[1]; ++i)
    for (int j = 0; j < nbr_points[0]; ++j, ++k)
      LLtoXYZ(x_points[j], y_points[i], coordinates[k]);

  *point_id = user_input_data_add_points(grid, location, coordinates);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_points_curve2d ( int const grid_id,
                               int const *nbr_points,
                               int const located,
                               double const *x_points,
                               double const *y_points,
                               int *point_id ) {

  enum yac_location location = yac_get_location(located);
  struct yac_basic_grid * grid = yac_unique_id_to_pointer(grid_id, "grid_id");

  size_t num_points = yac_basic_grid_get_data_size(grid, location);
  YAC_ASSERT_F(
    ((size_t)nbr_points[0] * (size_t)nbr_points[1]) == num_points,
    "ERROR(yac_cdef_points_curve2d): nbr_points does not match with grid"
    " (%zu*%zu != %zu)",
    (size_t)nbr_points[0], (size_t)nbr_points[1], num_points)

  check_x_vertices(
    x_points, (size_t)(nbr_points[0]), "yac_cdef_points_curve2d");
  check_y_vertices(
    y_points, (size_t)(nbr_points[1]), "yac_cdef_points_curve2d");

  coordinate_pointer coordinates =
    xmalloc(num_points * sizeof(*coordinates));
  for (int i = 0, k = 0; i < nbr_points[1]; ++i)
    for (int j = 0; j < nbr_points[0]; ++j, ++k)
      LLtoXYZ(x_points[k], y_points[k], coordinates[k]);

  *point_id = user_input_data_add_points(grid, location, coordinates);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_points_unstruct ( int const grid_id,
                                int const nbr_points,
                                int const located,
                                double const *x_points,
                                double const *y_points,
                                int *point_id ) {

  enum yac_location location = yac_get_location(located);
  struct yac_basic_grid * grid = yac_unique_id_to_pointer(grid_id, "grid_id");

  YAC_ASSERT(
    (size_t)nbr_points == yac_basic_grid_get_data_size(grid, location),
    "ERROR(yac_cdef_points_unstruct): nbr_points does not match with grid")

  check_x_vertices(x_points, (size_t)nbr_points, "yac_cdef_points_unstruct");
  check_y_vertices(y_points, (size_t)nbr_points, "yac_cdef_points_unstruct");

  coordinate_pointer coordinates =
    xmalloc((size_t)nbr_points * sizeof(*coordinates));
  for (int i = 0; i < nbr_points; ++i)
    LLtoXYZ(x_points[i], y_points[i], coordinates[i]);

  *point_id = user_input_data_add_points(grid, location, coordinates);
}

/* ---------------------------------------------------------------------- */

static int user_input_data_add_mask(
  struct yac_basic_grid * grid, enum yac_location location,
  int const * is_valid, size_t nbr_points, char const * name) {

  struct user_input_data_masks * curr_mask = xmalloc(1 * sizeof(*curr_mask));
  curr_mask->grid = grid;
  curr_mask->location = location;
  curr_mask->masks_idx =
    yac_basic_grid_add_mask(grid, location, is_valid, nbr_points, name);

  masks = xrealloc(masks, (num_masks + 1) * sizeof(*masks));
  masks[num_masks] = curr_mask;
  num_masks++;

  return yac_pointer_to_unique_id(curr_mask);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_mask_named ( int const grid_id,
                           int const nbr_points,
                           int const located,
                           int const * is_valid,
                           char const * name,
                           int *mask_id ) {

  enum yac_location location = yac_get_location(located);
  struct yac_basic_grid * grid = yac_unique_id_to_pointer(grid_id, "grid_id");

  YAC_ASSERT(
    (size_t)nbr_points == yac_basic_grid_get_data_size(grid, location),
    "ERROR(yac_cdef_mask_named): nbr_points does not match with grid")

  *mask_id =
    user_input_data_add_mask(
      grid, location, is_valid, (size_t)nbr_points, name);
}

void yac_cdef_mask ( int const grid_id,
                     int const nbr_points,
                     int const located,
                     int const * is_valid,
                     int *mask_id ) {

  yac_cdef_mask_named(grid_id, nbr_points, located, is_valid, NULL, mask_id);
}

void yac_cset_mask ( int const * is_valid, int points_id ) {

  struct user_input_data_points * points =
    yac_unique_id_to_pointer(points_id, "points_id");

  YAC_ASSERT(
    points->default_mask_id == INT_MAX,
    "ERROR(yac_cset_mask): default mask has already been set before")

  struct yac_basic_grid * grid = points->grid;
  enum yac_location location = points->location;
  size_t nbr_points = yac_basic_grid_get_data_size(grid, location);

  points->default_mask_id =
    user_input_data_add_mask(grid, location, is_valid, nbr_points, NULL);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_field_mask ( char const * name,
                           int const comp_id,
                           int const * point_ids,
                           int const * mask_ids,
                           int const num_pointsets,
                           int collection_size,
                           const char* timestep,
                           int time_unit,
                           int * field_id ) {

  YAC_ASSERT(
    num_pointsets >= 1,
    "ERROR(yac_cdef_field_mask): invalid number of pointsets")

  YAC_ASSERT(
    point_ids != NULL,"ERROR(yac_cdef_field_mask): no point_ids provided")

  struct yac_basic_grid * grid =
    ((struct user_input_data_points *)
      yac_unique_id_to_pointer(point_ids[0], "point_ids[0]"))->grid;

  struct interp_field interp_fields_buf;
  struct interp_field * interp_fields =
    (num_pointsets == 1)?
      &interp_fields_buf:
      xmalloc((size_t)num_pointsets * sizeof(*interp_fields));

  for (int i = 0; i < num_pointsets; ++i) {

    struct user_input_data_points * points =
      yac_unique_id_to_pointer(point_ids[i], "point_ids[i]");

    YAC_ASSERT(
      grid == points->grid,
      "ERROR(yac_cdef_field_mask): grid of point_ids do not match")

    size_t masks_idx;
    if (mask_ids[i] != INT_MAX) {

      struct user_input_data_masks * mask =
        yac_unique_id_to_pointer(mask_ids[i], "mask_ids[i]");

      YAC_ASSERT(
        mask->grid == points->grid,
        "ERROR(yac_cdef_field_mask): "
        "grids of mask and points do not match")
      YAC_ASSERT(
        mask->location == points->location,
        "ERROR(yac_cdef_field_mask): "
        "location of mask and points do not match")
      masks_idx = mask->masks_idx;
    } else {
      masks_idx = SIZE_MAX;
    }

    interp_fields[i].location = points->location;
    interp_fields[i].coordinates_idx = points->coordinates_idx;
    interp_fields[i].masks_idx = masks_idx;
  }

  struct user_input_data_component * comp_info =
    get_user_input_data_component(comp_id, "yac_cdef_field_mask");

  *field_id =
    yac_pointer_to_unique_id(
      yac_instance_add_field(
        comp_info->instance, name, comp_info->component_index,
        grid, interp_fields, num_pointsets, collection_size,
        yac_time_to_ISO(timestep, (enum yac_time_unit_type)time_unit)));

  if (num_pointsets > 1) free(interp_fields);
}

void yac_cdef_field ( char const * name,
                      int const comp_id,
                      int const * point_ids,
                      int const num_pointsets,
                      int collection_size,
                      const char* timestep,
                      int time_unit,
                      int * field_id ) {

  YAC_ASSERT(
    num_pointsets >= 1,
    "ERROR(yac_cdef_field): invalid number of pointsets")

  YAC_ASSERT(
    point_ids != NULL, "ERROR(yac_cdef_field): no point_ids provided")

  int * mask_ids = xmalloc((size_t)num_pointsets * sizeof(*mask_ids));

  for (int i = 0; i < num_pointsets; ++i)
    mask_ids[i] =
      ((struct user_input_data_points *)
         yac_unique_id_to_pointer(point_ids[i], "point_ids[i]"))->
           default_mask_id;

  yac_cdef_field_mask (
    name, comp_id, point_ids, mask_ids, num_pointsets,
    collection_size, timestep, time_unit, field_id);

  free(mask_ids);
}

void yac_cenable_field_frac_mask_instance(
  int yac_instance_id, const char* comp_name, const char* grid_name,
  const char* field_name, double frac_mask_fallback_value) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  yac_couple_config_field_enable_frac_mask(
    couple_config, comp_name, grid_name, field_name,
    frac_mask_fallback_value);
}

void yac_cenable_field_frac_mask(
  const char* comp_name, const char* grid_name, const char* field_name,
  double frac_mask_fallback_value) {
  check_default_instance_id("yac_cenable_field_frac_mask");
  yac_cenable_field_frac_mask_instance(
    default_instance_id, comp_name, grid_name, field_name,
    frac_mask_fallback_value);
}

void yac_cdef_component_metadata_instance(int yac_instance_id,
  const char* comp_name, const char* metadata) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config = yac_instance_get_couple_config(instance);
  yac_couple_config_component_set_metadata(couple_config, comp_name, metadata);
}

void yac_cdef_component_metadata(const char* comp_name, const char* metadata) {
  check_default_instance_id("yac_cdef_component_metadata");
  yac_cdef_component_metadata_instance(default_instance_id, comp_name, metadata);
}

void yac_cdef_grid_metadata_instance(int yac_instance_id, const char* grid_name,
  const char* metadata) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config
    = yac_instance_get_couple_config(instance);
  yac_couple_config_grid_set_metadata(couple_config, grid_name, metadata);
}

void yac_cdef_grid_metadata(const char* grid_name, const char* metadata) {
  check_default_instance_id("yac_cdef_grid_metadata");
  yac_cdef_grid_metadata_instance(default_instance_id, grid_name, metadata);
}

void yac_cdef_field_metadata_instance(int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name, const char* metadata) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config
    = yac_instance_get_couple_config(instance);
  yac_couple_config_field_set_metadata(couple_config, comp_name, grid_name, field_name, metadata);
}

void yac_cdef_field_metadata(const char* comp_name, const char* grid_name,
  const char* field_name, const char* metadata) {
  check_default_instance_id("yac_cdef_field_metadata");
  yac_cdef_field_metadata_instance(default_instance_id, comp_name, grid_name,
    field_name, metadata);
}

const char* yac_cget_component_metadata_instance(int yac_instance_id,
  const char* comp_name) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config = yac_instance_get_couple_config(instance);
  return yac_couple_config_component_get_metadata(couple_config, comp_name);
}

const char* yac_cget_component_metadata(const char* comp_name) {
  check_default_instance_id("yac_cget_component_metadata");
  return yac_cget_component_metadata_instance(default_instance_id, comp_name);
}

const char* yac_cget_grid_metadata_instance(int yac_instance_id,
  const char* grid_name) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config = yac_instance_get_couple_config(instance);
  return yac_couple_config_grid_get_metadata(couple_config, grid_name);
}

const char* yac_cget_grid_metadata(const char* grid_name) {
  check_default_instance_id("yac_cget_grid_metadata");
  return yac_cget_grid_metadata_instance(default_instance_id, grid_name);
}

const char* yac_cget_field_metadata_instance(int yac_instance_id,
  const char* comp_name, const char* grid_name, const char* field_name) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config = yac_instance_get_couple_config(instance);
  return yac_couple_config_field_get_metadata(couple_config, comp_name, grid_name,
    field_name);
}

const char* yac_cget_field_metadata(const char* comp_name, const char* grid_name,
  const char* field_name) {
  check_default_instance_id("yac_cget_field_metadata");
  return yac_cget_field_metadata_instance(default_instance_id, comp_name,
    grid_name, field_name);
}

static void init_ext_couple_config(
  struct yac_ext_couple_config * ext_couple_config) {
  ext_couple_config->weight_file = NULL;
  ext_couple_config->mapping_side = 1;
  ext_couple_config->scale_factor = 1.0;
  ext_couple_config->scale_summand = 0.0;
  ext_couple_config->num_src_mask_names = 0;
  ext_couple_config->src_mask_names = NULL;
  ext_couple_config->tgt_mask_name = NULL;
}

void yac_cget_ext_couple_config(int * ext_couple_config_id) {

  struct yac_ext_couple_config * ext_couple_config =
    xmalloc(1 * sizeof(*ext_couple_config));
  init_ext_couple_config(ext_couple_config);
  *ext_couple_config_id = yac_pointer_to_unique_id(ext_couple_config);
}

void yac_cfree_ext_couple_config(int ext_couple_config_id) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  free(ext_couple_config->weight_file);
  for (size_t i = 0; i < ext_couple_config->num_src_mask_names; ++i)
    free(ext_couple_config->src_mask_names[i]);
  free(ext_couple_config->src_mask_names);
  free(ext_couple_config->tgt_mask_name);
  free(ext_couple_config);
}

void yac_cset_ext_couple_config_weight_file_(
  struct yac_ext_couple_config * ext_couple_config,
  char const * weight_file) {
  free(ext_couple_config->weight_file);
  ext_couple_config->weight_file =
    (weight_file != NULL)?strdup(weight_file):NULL;
}

void yac_cset_ext_couple_config_weight_file(int ext_couple_config_id,
                                            char const * weight_file) {
  yac_cset_ext_couple_config_weight_file_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    weight_file);
}

void yac_cget_ext_couple_config_weight_file(int ext_couple_config_id,
                                            char const ** weight_file) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *weight_file = ext_couple_config->weight_file;
}

void yac_cset_ext_couple_config_mapping_side_(
  struct yac_ext_couple_config * ext_couple_config,
  int mapping_side) {
  YAC_ASSERT_F(
    (mapping_side == 0) ||
    (mapping_side == 1),
    "ERROR(yac_cset_ext_couple_config_mapping_side_): "
    "\"%d\" is not a valid mapping side (has to be 0 or 1)",
    mapping_side);
  ext_couple_config->mapping_side = mapping_side;
}

void yac_cset_ext_couple_config_mapping_side(int ext_couple_config_id,
                                             int mapping_side) {
  yac_cset_ext_couple_config_mapping_side_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    mapping_side);
}

void yac_cget_ext_couple_config_mapping_side(int ext_couple_config_id,
                                             int * mapping_side) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *mapping_side = ext_couple_config->mapping_side;
}

void yac_cset_ext_couple_config_scale_factor_(
  struct yac_ext_couple_config * ext_couple_config,
  double scale_factor) {
  YAC_ASSERT_F(isnormal(scale_factor),
    "ERROR(yac_cset_ext_couple_config_scale_factor_): "
    "\"%lf\" is not a valid scale factor", scale_factor);
  ext_couple_config->scale_factor = scale_factor;
}

void yac_cset_ext_couple_config_scale_factor(int ext_couple_config_id,
                                             double scale_factor) {
  yac_cset_ext_couple_config_scale_factor_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    scale_factor);
}

void yac_cget_ext_couple_config_scale_factor(int ext_couple_config_id,
                                             double * scale_factor) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *scale_factor = ext_couple_config->scale_factor;
}

void yac_cset_ext_couple_config_scale_summand_(
  struct yac_ext_couple_config * ext_couple_config,
  double scale_summand) {
  YAC_ASSERT_F((scale_summand == 0.0) || isnormal(scale_summand),
    "ERROR(yac_cset_ext_couple_config_scale_summand_): "
    "\"%lf\" is not a valid scale summand", scale_summand);
  ext_couple_config->scale_summand = scale_summand;
}

void yac_cset_ext_couple_config_scale_summand(int ext_couple_config_id,
                                             double scale_summand) {
  yac_cset_ext_couple_config_scale_summand_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    scale_summand);
}

void yac_cget_ext_couple_config_scale_summand(int ext_couple_config_id,
                                             double * scale_summand) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *scale_summand = ext_couple_config->scale_summand;
}

void yac_cset_ext_couple_config_src_mask_names_(
  struct yac_ext_couple_config * ext_couple_config,
  size_t num_src_mask_names, char const * const * src_mask_names) {
  YAC_ASSERT_F(num_src_mask_names >= 0,
    "ERROR(yac_cset_ext_couple_config_src_mask_names_): "
    "\"%zu\" is not a valid number of source mask names",
    num_src_mask_names)
  if (ext_couple_config->num_src_mask_names > 0)
    for (size_t i = 0; i < ext_couple_config->num_src_mask_names; ++i)
      free(ext_couple_config->src_mask_names[i]);
  free(ext_couple_config->src_mask_names);
  ext_couple_config->num_src_mask_names = num_src_mask_names;
  if (num_src_mask_names > 0) {
    ext_couple_config->src_mask_names =
      xmalloc(num_src_mask_names * sizeof(*src_mask_names));
    for (size_t i = 0; i < num_src_mask_names; ++i)
      ext_couple_config->src_mask_names[i] = strdup(src_mask_names[i]);
  } else ext_couple_config->src_mask_names = NULL;
}

void yac_cset_ext_couple_config_src_mask_names(
  int ext_couple_config_id, size_t num_src_mask_names,
  char const * const * src_mask_names) {
  yac_cset_ext_couple_config_src_mask_names_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    num_src_mask_names, src_mask_names);
}

void yac_cget_ext_couple_config_src_mask_names(
  int ext_couple_config_id, size_t * num_src_mask_names,
  char const * const ** src_mask_names) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *num_src_mask_names = ext_couple_config->num_src_mask_names;
  *src_mask_names = (char const * const *)ext_couple_config->src_mask_names;
}

void yac_cset_ext_couple_config_tgt_mask_name_(
  struct yac_ext_couple_config * ext_couple_config,
  char const * tgt_mask_name) {
  free(ext_couple_config->tgt_mask_name);
  ext_couple_config->tgt_mask_name =
    (tgt_mask_name != NULL)?strdup(tgt_mask_name):NULL;
}

void yac_cset_ext_couple_config_tgt_mask_name(
  int ext_couple_config_id, char const * tgt_mask_name) {
  yac_cset_ext_couple_config_tgt_mask_name_(
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"),
    tgt_mask_name);
}

void yac_cget_ext_couple_config_mask_name(
  int ext_couple_config_id, char const ** tgt_mask_name) {
  struct yac_ext_couple_config * ext_couple_config =
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id");
  *tgt_mask_name = ext_couple_config->tgt_mask_name;
}

void yac_cdef_couple_custom_instance_(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  struct yac_ext_couple_config * ext_couple_config) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");
  char const * coupling_timestep_iso8601 =
    yac_time_to_ISO(coupling_timestep, (enum yac_time_unit_type)time_unit);
  yac_instance_def_couple(instance,
    src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name,
    coupling_timestep_iso8601, time_reduction,
    interp_stack_config, src_lag, tgt_lag,
    ext_couple_config->weight_file,
    ext_couple_config->mapping_side,
    ext_couple_config->scale_factor,
    ext_couple_config->scale_summand,
    ext_couple_config->num_src_mask_names,
    (char const * const *)ext_couple_config->src_mask_names,
    ext_couple_config->tgt_mask_name);
}

void yac_cdef_couple_custom_instance(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  int ext_couple_config_id) {

  yac_cdef_couple_custom_instance_(
    yac_instance_id, src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_timestep,
    time_unit, time_reduction, interp_stack_config_id, src_lag, tgt_lag,
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"));
}

void yac_cdef_couple_custom(
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  int ext_couple_config_id) {

  check_default_instance_id("yac_cdef_couple_custom");
  yac_cdef_couple_custom_instance_(default_instance_id,
    src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_timestep,
    time_unit, time_reduction, interp_stack_config_id, src_lag, tgt_lag,
    yac_unique_id_to_pointer(ext_couple_config_id, "ext_couple_config_id"));
}

void yac_cdef_couple_instance(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag) {

  struct yac_ext_couple_config ext_couple_config;
  init_ext_couple_config(&ext_couple_config);

  yac_cdef_couple_custom_instance_(
    yac_instance_id, src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_timestep,
    time_unit, time_reduction, interp_stack_config_id, src_lag, tgt_lag,
    &ext_couple_config);
}

void yac_cdef_couple(
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
    char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
    char const * coupling_timestep, int time_unit, int time_reduction,
    int interp_stack_config_id, int src_lag, int tgt_lag){

  check_default_instance_id("yac_cdef_couple");
  yac_cdef_couple_instance(default_instance_id,
    src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name,
    coupling_timestep, time_unit, time_reduction,  interp_stack_config_id,
    src_lag, tgt_lag);
}

void yac_cdef_couple_instance_(
  int yac_instance_id,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_timestep, int time_unit, int time_reduction,
  int interp_stack_config_id, int src_lag, int tgt_lag,
  char const * weight_file, int mapping_side,
  double scale_factor, double scale_summand,
  int num_src_mask_names, char const * const * src_mask_names,
  char const * tgt_mask_name) {

  YAC_ASSERT_F(
    num_src_mask_names >= 0,
    "ERROR(yac_cdef_couple_instance_): "
    "\"%d\" is not a valid number of source mask names", num_src_mask_names)

  struct yac_ext_couple_config ext_couple_config;
  init_ext_couple_config(&ext_couple_config);
  yac_cset_ext_couple_config_weight_file_(
    &ext_couple_config, weight_file);
  yac_cset_ext_couple_config_mapping_side_(
    &ext_couple_config, mapping_side);
  yac_cset_ext_couple_config_scale_factor_(
    &ext_couple_config, scale_factor);
  yac_cset_ext_couple_config_scale_summand_(
    &ext_couple_config, scale_summand);
  yac_cset_ext_couple_config_src_mask_names_(
    &ext_couple_config, (size_t)num_src_mask_names, src_mask_names);
  yac_cset_ext_couple_config_tgt_mask_name_(
    &ext_couple_config, tgt_mask_name);

  yac_cdef_couple_custom_instance_(
    yac_instance_id, src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_timestep,
    time_unit, time_reduction, interp_stack_config_id, src_lag, tgt_lag,
    &ext_couple_config);
}

void yac_cdef_couple_(
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
    char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
    char const * coupling_timestep, int time_unit, int time_reduction,
    int interp_stack_config_id, int src_lag, int tgt_lag,
    char const * weight_file, int mapping_side,
    double scale_factor, double scale_summand,
    int num_src_mask_names, char const * const * src_mask_names,
    char const * tgt_mask_name) {

  check_default_instance_id("yac_cdef_couple_");
  yac_cdef_couple_instance_(default_instance_id,
    src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name,
    coupling_timestep, time_unit, time_reduction,  interp_stack_config_id,
    src_lag, tgt_lag, weight_file, mapping_side,
    scale_factor, scale_summand,
    num_src_mask_names, src_mask_names, tgt_mask_name);
}

/* ---------------------------------------------------------------------- */

void yac_ccheck_field_dimensions ( int field_id,
                                   int collection_size,
                                   int num_interp_fields,
                                   int const * interp_field_sizes ) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  YAC_ASSERT_F(
    collection_size ==
    (int)yac_get_coupling_field_collection_size(cpl_field),
    "ERROR(yac_ccheck_field_dimensions): mismatching collection sizes "
    "for component %s grid %s field %s (%d != %d)",
    yac_component_get_name(
      yac_get_coupling_field_component(cpl_field)),
    yac_basic_grid_get_name(
      yac_coupling_field_get_basic_grid(cpl_field)),
    yac_get_coupling_field_name(cpl_field),
    collection_size, (int)yac_get_coupling_field_collection_size(cpl_field));

  if (num_interp_fields != -1) {

    YAC_ASSERT_F(
      (size_t)num_interp_fields ==
      yac_coupling_field_get_num_interp_fields(cpl_field),
    "ERROR(yac_ccheck_field_dimensions): mismatching number of interp fields "
    "for component %s grid %s field %s (%d != %zu)",
    yac_component_get_name(
      yac_get_coupling_field_component(cpl_field)),
    yac_basic_grid_get_name(
      yac_coupling_field_get_basic_grid(cpl_field)),
    yac_get_coupling_field_name(cpl_field),
    num_interp_fields, yac_coupling_field_get_num_interp_fields(cpl_field));

    if (interp_field_sizes) {
      for (int interp_field_idx = 0; interp_field_idx < num_interp_fields;
          ++interp_field_idx) {

        YAC_ASSERT_F(
          (size_t)interp_field_sizes[interp_field_idx] ==
          yac_coupling_field_get_data_size(
            cpl_field,
            yac_get_coupling_field_get_interp_field_location(
              cpl_field, (size_t)interp_field_idx)),
          "ERROR(yac_ccheck_field_dimensions): mismatching interp field size "
          "for component %s grid %s field %s interp_field_idx %d (%d != %zu)",
          yac_component_get_name(
            yac_get_coupling_field_component(cpl_field)),
          yac_basic_grid_get_name(
            yac_coupling_field_get_basic_grid(cpl_field)),
          yac_get_coupling_field_name(cpl_field),
          interp_field_idx,
          interp_field_sizes[interp_field_idx],
          yac_coupling_field_get_data_size(
            cpl_field,
            yac_get_coupling_field_get_interp_field_location(
              cpl_field, (size_t)interp_field_idx)));
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void yac_cget_action(int field_id, int * action) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");
  enum yac_field_exchange_type exchange_type =
    yac_get_coupling_field_exchange_type(cpl_field);

  YAC_ASSERT(
    (exchange_type == NOTHING) ||
    (exchange_type == SOURCE) ||
    (exchange_type == TARGET),
    "ERROR(yac_cget_action): invalid field exchange type")

  switch(exchange_type) {
    default:
    case(NOTHING): {

      *action =  YAC_ACTION_NONE;
      break;
    }
    case(TARGET): {

      enum yac_action_type event_action =
        yac_event_check(yac_get_coupling_field_get_op_event(cpl_field));
      *action = (int)((event_action == RESTART)?GET_FOR_RESTART:event_action);
      break;
    }
    case(SOURCE): {

      enum yac_action_type event_action = NONE;
      unsigned num_puts = yac_get_coupling_field_num_puts(cpl_field);
      for (unsigned put_idx = 0; put_idx < num_puts; ++put_idx)
        event_action =
          MAX(
            event_action,
            yac_event_check(
              yac_get_coupling_field_put_op_event(cpl_field, put_idx)));
      *action = (int)((event_action == RESTART)?PUT_FOR_RESTART:event_action);
      break;
    }
  }
}

const char* yac_cget_field_datetime(int field_id){
  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");
  return yac_coupling_field_get_datetime(cpl_field);
}

/* ---------------------------------------------------------------------- */

static void inline yac_cupdate_(
  struct coupling_field * cpl_field, struct event * event, int is_source) {

  YAC_ASSERT_F(
    yac_event_check(event) == NONE,
    "ERROR(yac_cupdate): current action of %s field \"%s\" of "
    "componente \"%s\" is not YAC_ACTION_NONE",
    ((is_source)?"source":"target"), yac_get_coupling_field_name(cpl_field),
    yac_component_get_name(yac_get_coupling_field_component(cpl_field)))
  yac_event_update(event);
}

void yac_cupdate(int field_id) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");
  enum yac_field_exchange_type exchange_type =
    yac_get_coupling_field_exchange_type(cpl_field);

  YAC_ASSERT(
    (exchange_type == NOTHING) ||
    (exchange_type == SOURCE) ||
    (exchange_type == TARGET),
    "ERROR(yac_cupdate): invalid field exchange type")

  switch(exchange_type) {
    default:
    case(NOTHING): {
      break;
    }
    case(TARGET): {

      yac_cupdate_(
        cpl_field, yac_get_coupling_field_get_op_event(cpl_field), 0);
      break;
    }
    case(SOURCE): {

      unsigned num_puts = yac_get_coupling_field_num_puts(cpl_field);
      for (unsigned put_idx = 0; put_idx < num_puts; ++put_idx)
        yac_cupdate_(
          cpl_field,
          yac_get_coupling_field_put_op_event(cpl_field, put_idx), 1);
      break;
    }
  }
}

/* ---------------------------------------------------------------------- */

static double ** get_recv_field_pointers(
  int const field_id,
  int const collection_size,
  double *recv_field) { // recv_field[collection_size][Points]

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  YAC_ASSERT(
    yac_coupling_field_get_num_interp_fields(cpl_field) == 1,
    "ERROR(get_recv_field_pointers): invalid number of interpolation fields "
    "(should be one)")

  enum yac_location location =
    yac_coupling_field_get_interp_fields(cpl_field)->location;
  size_t num_points =
    yac_coupling_field_get_data_size(cpl_field, location);

  double ** recv_field_ = xmalloc(collection_size * sizeof(*recv_field_));

  for (int i = 0; i < collection_size; ++i) {
    recv_field_[i] = recv_field;
    recv_field += num_points;
  }

  return recv_field_;
}

void yac_cget_ ( int const field_id,
                 int const collection_size,
                 double *recv_field,   // recv_field[collection_size][Points]
                 int    *info,
                 int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, 1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double ** recv_field_ =
    get_recv_field_pointers(field_id, collection_size, recv_field);

  yac_cget ( field_id, collection_size, recv_field_, info, ierr );

  free(recv_field_);
}

/* ---------------------------------------------------------------------- */

static struct interpolation * yac_cget_pre_processing(
  int const field_id, int collection_size, int *info, int *ierr) {

  *info = (int)NONE;
  *ierr = 0;

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  if (yac_get_coupling_field_exchange_type(cpl_field) != TARGET)
    return NULL;

  /* --------------------------------------------------------------------
     Check for restart and coupling events
     -------------------------------------------------------------------- */

  struct event * event = yac_get_coupling_field_get_op_event(cpl_field);
  enum yac_action_type action = yac_event_check(event);
  *info = (int)((action == RESTART)?GET_FOR_RESTART:action);

  /* --------------------------------------------------------------------
     add one model time step (as provided in coupling configuration) to
     the current event date
     -------------------------------------------------------------------- */

  yac_event_update ( event );

  /* ------------------------------------------------------------------
     return in case we are already beyond the end of the run
     ------------------------------------------------------------------ */

  if ( action == OUT_OF_BOUND ) {

      fputs("WARNING: yac get action is beyond end of run date!\n", stderr);
      fputs("WARNING: yac get action is beyond end of run date!\n", stdout);

      return  NULL;
  }

  /* --------------------------------------------------------------------
     start actions
     -------------------------------------------------------------------- */

  if ( action == NONE ) return NULL;

  YAC_ASSERT(
    collection_size == yac_get_coupling_field_collection_size(cpl_field),
    "ERROR: collection size does not match with coupling configuration.")

  return yac_get_coupling_field_get_op_interpolation(cpl_field);
}

void yac_cget ( int const field_id,
                int collection_size,
                double ** recv_field, // recv_field[collection_size][Points]
                int    *info,
                int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, 1, NULL);

  struct interpolation * interpolation =
    yac_cget_pre_processing(field_id, collection_size, info, ierr);

  if ((*ierr == 0) && (interpolation != NULL))
    yac_interpolation_execute_get(interpolation, recv_field);
}

/* ---------------------------------------------------------------------- */

static void get_send_field_pointers(
  int const field_id,
  int const collection_size,
  double * send_field, // send_field[collection_size][nPointSets][Points]
  double * send_frac_mask,
  double **** send_field_,
  double **** send_frac_mask_) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  size_t num_interp_fields =
    yac_coupling_field_get_num_interp_fields(cpl_field);

  size_t * num_points = xmalloc(num_interp_fields * sizeof(*num_points));

  for (size_t i = 0; i < num_interp_fields; i++) {

    enum yac_location location =
      yac_coupling_field_get_interp_fields(cpl_field)[i].location;
    num_points[i] = yac_coupling_field_get_data_size(cpl_field, location);
  }

  *send_field_ = xmalloc(collection_size * sizeof(**send_field_));
  for (int i = 0; i < collection_size; ++i) {
    (*send_field_)[i] = xmalloc(num_interp_fields * sizeof(***send_field_));
    for (size_t j = 0; j < num_interp_fields; ++j) {
      (*send_field_)[i][j] = send_field;
      send_field += num_points[j];
    }
  }
  if (send_frac_mask != NULL) {
    *send_frac_mask_ =
      xmalloc(collection_size * sizeof(**send_frac_mask_));
    for (int i = 0; i < collection_size; ++i) {
      (*send_frac_mask_)[i] =
        xmalloc(num_interp_fields * sizeof(***send_frac_mask_));
      for (size_t j = 0; j < num_interp_fields; ++j) {
        (*send_frac_mask_)[i][j] = send_frac_mask;
        send_frac_mask += num_points[j];
      }
    }
  }
  free(num_points);
}

void yac_cput_ ( int const field_id,
                 int const collection_size,
                 double *send_field,   // send_field[collection_size]
                                       //           [nPointSets]
                                       //           [Points]
                 int    *info,
                 int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  get_send_field_pointers(
    field_id, collection_size, send_field, NULL, &send_field_, NULL);

  yac_cput ( field_id, collection_size, send_field_, info, ierr );

  for (int i = 0; i < collection_size; ++i) free(send_field_[i]);
  free(send_field_);
}

void yac_cput_frac_ ( int const field_id,
                      int const collection_size,
                      double *send_field,   // send_field[collection_size]
                                            //           [nPointSets]
                                            //           [Points]
                      double *send_frac_mask,   // send_frac_mask[collection_size]
                                                //               [nPointSets]
                                                //               [Points]
                      int    *info,
                      int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  double *** send_frac_mask_;
  get_send_field_pointers(
    field_id, collection_size, send_field, send_frac_mask,
    &send_field_, &send_frac_mask_);

  yac_cput_frac(
    field_id, collection_size, send_field_, send_frac_mask_, info, ierr);

  for (int i = 0; i < collection_size; ++i) {
    free(send_field_[i]);
    free(send_frac_mask_[i]);
  }
  free(send_field_);
  free(send_frac_mask_);
}

/* ---------------------------------------------------------------------- */

static void get_send_field_pointers_ptr_(
  int const field_id,
  int const collection_size,
  double **send_field, // send_field[collection_size*nPointSets][Points]
  double **send_frac_mask,
  double ****send_field_,
  double ****send_frac_mask_) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  size_t num_interp_fields =
    yac_coupling_field_get_num_interp_fields(cpl_field);

  size_t * num_points = xmalloc(num_interp_fields * sizeof(*num_points));

  for (size_t i = 0; i < num_interp_fields; i++) {

    enum yac_location location =
      yac_coupling_field_get_interp_fields(cpl_field)[i].location;
    num_points[i] = yac_coupling_field_get_data_size(cpl_field, location);
  }

  *send_field_ = xmalloc(collection_size * sizeof(**send_field_));
  for (int i = 0, k = 0; i < collection_size; ++i) {

    (*send_field_)[i] = xmalloc(num_interp_fields * sizeof(***send_field_));
    for (size_t j = 0; j < num_interp_fields; ++j, ++k)
      (*send_field_)[i][j] = send_field[k];
  }
  if (send_frac_mask != NULL) {
    *send_frac_mask_ = xmalloc(collection_size * sizeof(**send_frac_mask_));
    for (int i = 0, k = 0; i < collection_size; ++i) {

      (*send_frac_mask_)[i] =
        xmalloc(num_interp_fields * sizeof(***send_frac_mask_));
      for (size_t j = 0; j < num_interp_fields; ++j, ++k)
        (*send_frac_mask_)[i][j] = send_frac_mask[k];
    }
  }

  free(num_points);
}

void yac_cput_ptr_ ( int const field_id,
                     int const collection_size,
                     double ** send_field,   // send_field[collection_size *
                                             //            nPointSets][Points]
                     int    *info,
                     int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  get_send_field_pointers_ptr_(
    field_id, collection_size, send_field, NULL, &send_field_, NULL);

  yac_cput ( field_id, collection_size, send_field_, info, ierr );

  for (int i = 0; i < collection_size; ++i) free(send_field_[i]);
  free(send_field_);
}

void yac_cput_frac_ptr_ ( int const field_id,
                          int const collection_size,
                          double ** send_field,   // send_field[collection_size *
                                                  //            nPointSets][Points]
                          double ** send_frac_mask, // send_frac_mask[collection_size *
                                                    //                nPointSets][Points]
                          int    *info,
                          int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  double *** send_frac_mask_;
  get_send_field_pointers_ptr_(
    field_id, collection_size, send_field, send_frac_mask,
    &send_field_, &send_frac_mask_);

  yac_cput_frac(
    field_id, collection_size, send_field_, send_frac_mask_, info, ierr);

  for (int i = 0; i < collection_size; ++i) {
    free(send_field_[i]);
    free(send_frac_mask_[i]);
  }
  free(send_field_);
  free(send_frac_mask_);
}

/* ---------------------------------------------------------------------- */

void yac_ctest(int field_id, int * flag) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  if (yac_get_coupling_field_exchange_type(cpl_field) != SOURCE) {
    *flag = 1;
    return;
  }

  *flag = 0;
  for (
    unsigned put_idx = 0;
    (put_idx < yac_get_coupling_field_num_puts(cpl_field)) && !*flag;
    ++put_idx)
    *flag =
      yac_interpolation_execute_test(
        yac_get_coupling_field_put_op_interpolation(cpl_field, put_idx));
}

/* ---------------------------------------------------------------------- */

void * yac_get_field_put_mask_c2f(int field_id) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  return (void*)yac_get_coupling_field_put_mask(cpl_field);
}

void * yac_get_field_get_mask_c2f(int field_id) {

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  return (void*)yac_get_coupling_field_get_mask(cpl_field);
}

/* ---------------------------------------------------------------------- */

static struct interpolation * yac_cput_pre_processing(
  struct coupling_field * cpl_field, unsigned put_idx,
  int const collection_size, double *** send_field,
  double **** send_field_acc_, double *** send_frac_mask,
  double **** send_frac_mask_acc_, int *info, int *ierr) {

  *ierr = 0;

  int with_frac_mask = send_frac_mask != NULL;

  /* ------------------------------------------------------------------
     Check for restart and coupling events
     ------------------------------------------------------------------ */

  struct event * event = yac_get_coupling_field_put_op_event(cpl_field, put_idx);
  enum yac_action_type action = yac_event_check(event);
  *info = (int)((action == RESTART)?PUT_FOR_RESTART:action);

  /* ------------------------------------------------------------------
     add one model time step (as provided in coupling configuration) to
     the current event date
     ------------------------------------------------------------------ */

  yac_event_update(event);

  /* ------------------------------------------------------------------
     return in case we are already beyond the end of the run
     ------------------------------------------------------------------ */

  if (action == OUT_OF_BOUND) {

    fputs("WARNING: yac put action is beyond end of run date!\n", stderr);
    fputs("WARNING: yac put action is beyond end of run date!\n", stdout);

    *send_field_acc_ = NULL;
    *send_frac_mask_acc_ = NULL;
    return NULL;
  }

  /* ------------------------------------------------------------------
     start actions
     ------------------------------------------------------------------ */

  if ( action == NONE ) return NULL;

  /* ------------------------------------------------------------------
     If it is time for restart, set appropriate flags and continue
     ------------------------------------------------------------------ */

  /* ------------------------------------------------------------------
     First deal with instant sends
     ------------------------------------------------------------------ */

  int ** put_mask = yac_get_coupling_field_put_mask(cpl_field);
  size_t num_interp_fields =
    yac_coupling_field_get_num_interp_fields(cpl_field);

  int time_operation = yac_get_event_time_operation(event);
  if ( time_operation == TIME_NONE ) {

    *info = MAX((int)COUPLING, *info);

    if (with_frac_mask) {

      // apply fractional mask to send field
      double *** send_field_acc =
        yac_get_coupling_field_put_op_send_field_acc(cpl_field, put_idx);
      *send_field_acc_ = send_field_acc;
      *send_frac_mask_acc_ = send_frac_mask;

      for (int h = 0; h < collection_size; h++) {
        for (size_t i = 0; i < num_interp_fields; i++) {
          size_t num_points =
            yac_coupling_field_get_data_size(
              cpl_field,
              yac_coupling_field_get_interp_fields(cpl_field)[i].location);
          if (put_mask) {
            for (size_t j = 0; j < num_points; ++j) {
              if (put_mask[i][j]) {
                double frac_mask = send_frac_mask[h][i][j];
                if (frac_mask != 0.0)
                  send_field_acc[h][i][j] = send_field[h][i][j] * frac_mask;
              }
            }
          } else {
            for (size_t j = 0; j < num_points; ++j) {
              double frac_mask = send_frac_mask[h][i][j];
              if (frac_mask != 0.0)
                send_field_acc[h][i][j] = send_field[h][i][j] * frac_mask;
            }
          }
        }
      }
    } else {
      *send_field_acc_ = NULL;
      *send_frac_mask_acc_ = NULL;
    }
    return yac_get_coupling_field_put_op_interpolation(cpl_field, put_idx);
  }

  /* ------------------------------------------------------------------
     Accumulation & Averaging
     ------------------------------------------------------------------ */

  YAC_ASSERT(
    (time_operation == TIME_ACCUMULATE) ||
    (time_operation == TIME_AVERAGE) ||
    (time_operation == TIME_MINIMUM) ||
    (time_operation == TIME_MAXIMUM),
    "ERROR(yac_cput_pre_processing): invalid time operation type")

  int time_accumulation_count =
    yac_get_coupling_field_put_op_time_accumulation_count(cpl_field, put_idx);

  time_accumulation_count++;

  // if this is the first accumulation step
  if (time_accumulation_count == 1) {

    double send_field_acc_init_value;
    switch (time_operation) {
      default:
      case(TIME_ACCUMULATE):
      case(TIME_AVERAGE):
        send_field_acc_init_value = 0.0;
        break;
      case(TIME_MINIMUM):
        send_field_acc_init_value = DBL_MAX;
        break;
      case(TIME_MAXIMUM):
        send_field_acc_init_value = -DBL_MAX;
        break;
    }

    /* initalise memory */

    yac_init_coupling_field_put_op_send_field_acc(
      cpl_field, put_idx, send_field_acc_init_value);
    if (with_frac_mask)
      yac_init_coupling_field_put_op_send_frac_mask_acc(
        cpl_field, put_idx, 0.0);
  }

  /* accumulate data */

  double *** send_field_acc =
    yac_get_coupling_field_put_op_send_field_acc(cpl_field, put_idx);
  *send_field_acc_ = send_field_acc;
  double *** send_frac_mask_acc;
  if (with_frac_mask) {
    send_frac_mask_acc =
      yac_get_coupling_field_put_op_send_frac_mask_acc(cpl_field, put_idx);
    *send_frac_mask_acc_ = send_frac_mask_acc;
  } else {
    send_frac_mask_acc = NULL;
    *send_frac_mask_acc_ = NULL;
  }

#define NO_CHECK (1)
#define PUT_CHECK (put_mask[i][j])
#define SUM +=
#define ASSIGN =
#define AGGREGATE_FRAC(CHECK, EXTRA_CHECK, ACCU_OP) \
  { \
    for (size_t j = 0; j < num_points; ++j) { \
      double frac_mask = send_frac_mask[h][i][j]; \
      double send_field_value = send_field[h][i][j] * frac_mask; \
      if (CHECK && (EXTRA_CHECK)) { \
        if (frac_mask != 0.0) { \
          send_field_acc[h][i][j] ACCU_OP send_field_value; \
          send_frac_mask_acc[h][i][j] ACCU_OP frac_mask; \
        } \
      } \
    } \
  }
#define AGGREATE_NOFRAC(CHECK, EXTRA_CHECK, ACCU_OP) \
  { \
    for (size_t j = 0; j < num_points; ++j) {\
      double send_field_value = send_field[h][i][j]; \
      if (CHECK && (EXTRA_CHECK)) \
        send_field_acc[h][i][j] ACCU_OP send_field_value; \
    } \
  }
#define AGGREGATE(EXTRA_CHECK, ACCU_OP) \
  { \
    for (int h = 0; h < collection_size; h++) { \
      for (size_t i = 0; i < num_interp_fields; i++) { \
        size_t num_points = \
          yac_coupling_field_get_data_size( \
            cpl_field, \
            yac_coupling_field_get_interp_fields(cpl_field)[i].location); \
        if (with_frac_mask) { \
          if (put_mask) AGGREGATE_FRAC(PUT_CHECK, EXTRA_CHECK, ACCU_OP) \
          else          AGGREGATE_FRAC(NO_CHECK,  EXTRA_CHECK, ACCU_OP) \
        } else { \
          if (put_mask) AGGREATE_NOFRAC(PUT_CHECK, EXTRA_CHECK, ACCU_OP) \
          else          AGGREATE_NOFRAC(NO_CHECK,  EXTRA_CHECK, ACCU_OP) \
        } \
      } \
    } \
    break; \
  }

  switch (time_operation) {
    default:
    case(TIME_ACCUMULATE):
    case(TIME_AVERAGE):
      AGGREGATE(NO_CHECK, SUM)
    case(TIME_MINIMUM):
      AGGREGATE(send_field_acc[h][i][j] > send_field_value, ASSIGN)
    case(TIME_MAXIMUM):
      AGGREGATE(send_field_acc[h][i][j] < send_field_value, ASSIGN)
  }

#undef AGGREGATE
#undef AGGREATE_NOFRAC
#undef AGGREGATE_FRAC
#undef ASSIGN
#undef SUM
#undef PUT_CHECK
#undef NO_CHECK

/* --------------------------------------------------------------------
   Check whether we have to perform the coupling in this call
   -------------------------------------------------------------------- */

  if (action == REDUCTION) {
    yac_set_coupling_field_put_op_time_accumulation_count(
      cpl_field, put_idx, time_accumulation_count);
    return NULL;
  }

/* --------------------------------------------------------------------
   Average data if required
   -------------------------------------------------------------------- */

  if (( time_operation == TIME_AVERAGE ) ||
      ( time_operation == TIME_ACCUMULATE )) {

    double weight = 1.0 / (double)time_accumulation_count;

#define NO_CHECK (1)
#define PUT_CHECK (put_mask[i][j])
#define WEIGHT_ACC_(ACC, CHECK, EXTRA_CHECK) \
  { \
    for (size_t j = 0; j < num_points; j++) \
      if (CHECK && (EXTRA_CHECK)) ACC[h][i][j] *= weight; \
  }
#define WEIGHT_ACC(ACC, EXTRA_CHECK) \
  { \
    if (put_mask) WEIGHT_ACC_(ACC, PUT_CHECK, EXTRA_CHECK) \
    else          WEIGHT_ACC_(ACC, NO_CHECK,  EXTRA_CHECK) \
  }

    for (int h = 0; h < collection_size; ++h) {
      for (size_t i = 0; i < num_interp_fields; i++) {
        size_t num_points =
          yac_coupling_field_get_data_size(
            cpl_field,
            yac_coupling_field_get_interp_fields(cpl_field)[i].location);
        if (with_frac_mask)
          WEIGHT_ACC(send_frac_mask_acc, NO_CHECK)
        if (time_operation == TIME_AVERAGE) {
          if (with_frac_mask)
            WEIGHT_ACC(send_field_acc, send_frac_mask_acc[h][i][j] != 0.0)
          else
            WEIGHT_ACC(send_field_acc, NO_CHECK)
        }
      }
    }
  }

#undef NO_CHECK
#undef PUT_CHECK
#undef WEIGHT_ACC_
#undef WEIGHT_ACC

  // reset time_accumulation_count
  yac_set_coupling_field_put_op_time_accumulation_count(
    cpl_field, put_idx, 0);

/* --------------------------------------------------------------------
   return interpolation
   -------------------------------------------------------------------- */

  *info = MAX((int)COUPLING, *info);
  return yac_get_coupling_field_put_op_interpolation(cpl_field, put_idx);
}

void yac_cput_frac ( int const field_id,
                     int const collection_size,
                     double *** const send_field,     // send_field[collection_size]
                                                      //           [nPointSets][Points]
                     double *** const send_frac_mask, // send_frac_mask[collection_size]
                                                      //               [nPointSets][Points]
                     int    *info,
                     int    *ierr ) {

  yac_ccheck_field_dimensions(field_id, collection_size, -1, NULL);

  *info = NONE;
  *ierr = 0;

  struct coupling_field * cpl_field =
    yac_unique_id_to_pointer(field_id, "field_id");

  if (yac_get_coupling_field_exchange_type(cpl_field) != SOURCE)
    return;

  YAC_ASSERT(
    collection_size == yac_get_coupling_field_collection_size(cpl_field),
    "ERROR: collection size does not match with coupling configuration.")

  for (unsigned put_idx = 0;
       put_idx < yac_get_coupling_field_num_puts(cpl_field); ++put_idx) {

    int curr_action;
    int curr_ierr;
    double *** send_field_acc;
    double *** send_frac_mask_acc;
    struct interpolation * interpolation =
      yac_cput_pre_processing(
        cpl_field, put_idx, collection_size,
        send_field, &send_field_acc,
        send_frac_mask, &send_frac_mask_acc,
        &curr_action, &curr_ierr);

    *info = MAX(*info, curr_action);
    *ierr = MAX(*ierr, curr_ierr);

    /* ------------------------------------------------------------------
       return in case we are already beyond the end of the run
       ------------------------------------------------------------------ */
    if (curr_action == OUT_OF_BOUND) {
      *info = OUT_OF_BOUND;
      return;
    }

    /* ------------------------------------------------------------------
       in case there is nothing to be done for the current put
       ------------------------------------------------------------------ */
    if ((curr_action == NONE) || (curr_action == REDUCTION)) continue;

    /* ------------------------------------------------------------------
       in case we are supposed to couple
       ------------------------------------------------------------------ */
    if ((*ierr == 0) && (interpolation != NULL)) {

      int with_frac_mask =
        yac_interpolation_with_frac_mask(interpolation);

      YAC_ASSERT_F(
        (with_frac_mask && (send_frac_mask != NULL)) ||
        !with_frac_mask,
        "ERROR: interpolation for field \"%s\" was built for dynamic "
        "fractional masking, but no mask was provided",
        yac_get_coupling_field_name(cpl_field));
      YAC_ASSERT_F(
        (!with_frac_mask && (send_frac_mask == NULL)) ||
        with_frac_mask,
        "ERROR: interpolation for field \"%s\" was not built for dynamic "
        "fractional masking, but a mask was provided",
        yac_get_coupling_field_name(cpl_field));

      double *** send_field_ptr =
        (send_field_acc == NULL)?send_field:send_field_acc;
      double *** send_frac_mask_ptr =
        (send_frac_mask_acc == NULL)?send_frac_mask:send_frac_mask_acc;

      if (with_frac_mask)
        yac_interpolation_execute_put_frac(
          interpolation, send_field_ptr, send_frac_mask_ptr);
      else
        yac_interpolation_execute_put(interpolation, send_field_ptr);
    }
  }
}

void yac_cput ( int const field_id,
                int const collection_size,
                double *** const send_field,   // send_field[collection_size]
                                               //           [nPointSets][Points]
                int    *info,
                int    *ierr ) {

  yac_cput_frac(
    field_id, collection_size, send_field, NULL, info, ierr);
}

/* ---------------------------------------------------------------------- */


void yac_cexchange_ ( int const send_field_id,
                      int const recv_field_id,
                      int const collection_size,
                      double *send_field, // send_field[collection_size]
                                          //           [nPointSets]
                                          //           [Points]
                      double *recv_field, // recv_field[collection_size][Points]
                      int    *send_info,
                      int    *recv_info,
                      int    *ierr ) {

  yac_ccheck_field_dimensions(send_field_id, collection_size, 1, NULL);
  yac_ccheck_field_dimensions(recv_field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  get_send_field_pointers(
    send_field_id, collection_size, send_field, NULL, &send_field_, NULL);
  double ** recv_field_ =
    get_recv_field_pointers(recv_field_id, collection_size, recv_field);

  yac_cexchange(
    send_field_id, recv_field_id, collection_size, send_field_, recv_field_,
    send_info, recv_info, ierr);

  free(recv_field_);
  for (int i = 0; i < collection_size; ++i) free(send_field_[i]);
  free(send_field_);
}

/* ---------------------------------------------------------------------- */


void yac_cexchange_frac_ ( int const send_field_id,
                           int const recv_field_id,
                           int const collection_size,
                           double *send_field, // send_field[collection_size]
                                               //           [nPointSets]
                                               //           [Points]
                           double *send_frac_mask, // send_frac_mask[collection_size]
                                                   //               [nPointSets]
                                                   //               [Points]
                           double *recv_field, // recv_field[collection_size][Points]
                           int    *send_info,
                           int    *recv_info,
                           int    *ierr ) {

  yac_ccheck_field_dimensions(send_field_id, collection_size, 1, NULL);
  yac_ccheck_field_dimensions(recv_field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  double *** send_frac_mask_;
  get_send_field_pointers(
    send_field_id, collection_size, send_field, send_frac_mask,
    &send_field_, &send_frac_mask_);
  double ** recv_field_ =
    get_recv_field_pointers(recv_field_id, collection_size, recv_field);

  yac_cexchange_frac(
    send_field_id, recv_field_id, collection_size, send_field_, send_frac_mask_,
    recv_field_, send_info, recv_info, ierr);

  free(recv_field_);
  for (int i = 0; i < collection_size; ++i) {
    free(send_field_[i]);
    free(send_frac_mask_[i]);
  }
  free(send_field_);
  free(send_frac_mask_);
}

/* ---------------------------------------------------------------------- */


void yac_cexchange_ptr_ ( int const send_field_id,
                          int const recv_field_id,
                          int const collection_size,
                          double ** send_field, // send_field[collection_size *
                                                //            nPointSets][Points]
                          double ** recv_field, // recv_field[collection_size][Points]
                          int    *send_info,
                          int    *recv_info,
                          int    *ierr ) {

  yac_ccheck_field_dimensions(send_field_id, collection_size, 1, NULL);
  yac_ccheck_field_dimensions(recv_field_id, collection_size, -1, NULL);


  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  get_send_field_pointers_ptr_(
    send_field_id, collection_size, send_field, NULL, &send_field_, NULL);

  yac_cexchange(
    send_field_id, recv_field_id, collection_size, send_field_, recv_field,
    send_info, recv_info, ierr);

  for (int i = 0; i < collection_size; ++i) free(send_field_[i]);
  free(send_field_);
}

/* ---------------------------------------------------------------------- */


void yac_cexchange_frac_ptr_ ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double ** send_field, // send_field[collection_size *
                                                     //            nPointSets][Points]
                               double ** send_frac_mask, // send_frac_mask[collection_size *
                                                         //                nPointSets][Points]
                               double ** recv_field, // recv_field[collection_size][Points]
                               int    *send_info,
                               int    *recv_info,
                               int    *ierr ) {

  yac_ccheck_field_dimensions(send_field_id, collection_size, 1, NULL);
  yac_ccheck_field_dimensions(recv_field_id, collection_size, -1, NULL);

  /* Needed to transfer from Fortran data structure to C */
  double *** send_field_;
  double *** send_frac_mask_;
  get_send_field_pointers_ptr_(
    send_field_id, collection_size, send_field, send_frac_mask,
    &send_field_, &send_frac_mask_);

  yac_cexchange_frac(
    send_field_id, recv_field_id, collection_size, send_field_, send_frac_mask_,
    recv_field, send_info, recv_info, ierr);

  for (int i = 0; i < collection_size; ++i) {
    free(send_field_[i]);
    free(send_frac_mask_[i]);
  }
  free(send_field_);
  free(send_frac_mask_);
}

/* ---------------------------------------------------------------------- */

void yac_cexchange_frac ( int const send_field_id,
                          int const recv_field_id,
                          int const collection_size,
                          double *** const send_field, // send_field[collection_size]
                                                       //           [nPointSets][Points]
                          double *** const send_frac_mask, // send_frac_mask[collection_size]
                                                           //               [nPointSets][Points]
                          double ** recv_field, // recv_field[collection_size][Points]
                          int    *send_info,
                          int    *recv_info,
                          int    *ierr ) {

  yac_ccheck_field_dimensions(send_field_id, collection_size, 1, NULL);
  yac_ccheck_field_dimensions(recv_field_id, collection_size, -1, NULL);

  *send_info = NONE;
  *ierr = -1;

  struct coupling_field * send_cpl_field =
    yac_unique_id_to_pointer(send_field_id, "send_field_id");

  if (yac_get_coupling_field_exchange_type(send_cpl_field) != SOURCE) {
    yac_cget(recv_field_id, collection_size, recv_field, recv_info, ierr);
    return;
  }

  YAC_ASSERT(
    collection_size == yac_get_coupling_field_collection_size(send_cpl_field),
    "ERROR(yac_cexchange_frac): "
    "collection size does not match with coupling configuration.")

  YAC_ASSERT(
    yac_get_coupling_field_num_puts(send_cpl_field) == 1,
    "ERROR(yac_cexchange_frac): more than one put per field is not supported "
    "for yac_cexchange_frac.")

  int send_action;
  int send_ierr;
  double *** send_field_acc;
  double *** send_frac_mask_acc;
  struct interpolation * put_interpolation =
    yac_cput_pre_processing(
      send_cpl_field, 0, collection_size, send_field,
      &send_field_acc, send_frac_mask, &send_frac_mask_acc,
      &send_action, &send_ierr);

  *send_info = MAX(*send_info, send_action);
  *ierr = MAX(*ierr, send_ierr);

  /* ------------------------------------------------------------------
     return in case we are already beyond the end of the run for the puts
     or there is nothing to be done for the current put
     ------------------------------------------------------------------ */
  if ((send_action == OUT_OF_BOUND) ||
      (send_action == NONE) ||
      (send_action == REDUCTION)) {
    yac_cget(recv_field_id, collection_size, recv_field, recv_info, ierr);
    return;
  }

  /* ------------------------------------------------------------------
     check the get
     ------------------------------------------------------------------ */

  struct interpolation * get_interpolation =
    yac_cget_pre_processing(recv_field_id, collection_size, recv_info, ierr);

  /* ------------------------------------------------------------------
     do the required exchanges
     ------------------------------------------------------------------ */

  if (*ierr == 0) {

    double *** send_field_ptr =
      (send_field_acc == NULL)?send_field:send_field_acc;
    double *** send_frac_mask_ptr =
      (send_frac_mask_acc == NULL)?send_frac_mask:send_frac_mask_acc;

    int with_frac_mask =
      yac_interpolation_with_frac_mask(put_interpolation);

    // if the get is active
    if (get_interpolation != NULL) {

      YAC_ASSERT(
        get_interpolation == put_interpolation,
        "ERROR(yac_cexchange): send_field_id and recv_field_id do not match")

      YAC_ASSERT_F(
        (with_frac_mask && (send_frac_mask != NULL)) ||
        !with_frac_mask,
        "ERROR: interpolation for field \"%s\" was built for dynamic "
        "fractional masking, but no mask was provided",
        yac_get_coupling_field_name(send_cpl_field));
      YAC_ASSERT_F(
        (!with_frac_mask && (send_frac_mask == NULL)) ||
        with_frac_mask,
        "ERROR: interpolation for field \"%s\" was not built for dynamic "
        "fractional masking, but a mask was provided",
        yac_get_coupling_field_name(send_cpl_field));

      if (with_frac_mask)
        yac_interpolation_execute_frac(
          put_interpolation, send_field_ptr, send_frac_mask_ptr, recv_field);
      else
        yac_interpolation_execute(
          put_interpolation, send_field_ptr, recv_field);

    } else {

      YAC_ASSERT_F(
        (with_frac_mask && (send_frac_mask != NULL)) ||
        !with_frac_mask,
        "ERROR: interpolation for field \"%s\" was built for dynamic "
        "fractional masking, but no mask was provided",
        yac_get_coupling_field_name(send_cpl_field));
      YAC_ASSERT_F(
        (!with_frac_mask && (send_frac_mask == NULL)) ||
        with_frac_mask,
        "ERROR: interpolation for field \"%s\" was not built for dynamic "
        "fractional masking, but a mask was provided",
        yac_get_coupling_field_name(send_cpl_field));

      // just execute the put
      if (with_frac_mask)
        yac_interpolation_execute_put_frac(
          put_interpolation, send_field_ptr, send_frac_mask_ptr);
      else
        yac_interpolation_execute_put(put_interpolation, send_field_ptr);
    }
  }
}



void yac_cexchange ( int const send_field_id,
                     int const recv_field_id,
                     int const collection_size,
                     double *** const send_field, // send_field[collection_size]
                                                  //           [nPointSets][Points]
                     double ** recv_field, // recv_field[collection_size][Points]
                     int    *send_info,
                     int    *recv_info,
                     int    *ierr ) {

  yac_cexchange_frac(
    send_field_id, recv_field_id, collection_size,
    send_field, NULL, recv_field, send_info, recv_info, ierr);
}

/* ---------------------------------------------------------------------- */

void yac_csync_def_instance ( int yac_instance_id ){
  yac_instance_sync_def(
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"));
}

void yac_csync_def ( void ){
  check_default_instance_id("yac_csync");
  yac_csync_def_instance(default_instance_id);
}


/* ---------------------------------------------------------------------- */

void yac_cenddef_instance(int yac_instance_id) {

  yac_instance_setup(
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"));
}

void yac_cenddef (void) {

  check_default_instance_id("yac_cenddef");
  yac_cenddef_instance(default_instance_id);
}

void yac_cenddef_and_emit_config_instance(
  int yac_instance_id, int emit_flags, char ** config) {

  *config =
    yac_instance_setup_and_emit_config(
      yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id"),
      emit_flags);
}

void yac_cenddef_and_emit_config ( int emit_flags, char ** config ) {

  check_default_instance_id("yac_cenddef");
  yac_cenddef_and_emit_config_instance(
    default_instance_id, emit_flags, config);
}

/* ----------------------------------------------------------------------
                   query functions
   ----------------------------------------------------------------------*/

int yac_cget_field_id_instance(int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name){
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct coupling_field* cpl_field =
    yac_instance_get_field(instance, comp_name, grid_name, field_name);
  YAC_ASSERT_F(cpl_field != NULL,
    "ERROR(yac_cget_field_id_instance): "
    "no field '%s' defined on the local process for "
    "component '%s' and grid '%s'", field_name, comp_name, grid_name);
  return yac_lookup_pointer(cpl_field);
}

int yac_cget_field_id(const char* comp_name, const char* grid_name, const char* field_name){
  check_default_instance_id("yac_cget_field_id");
  return yac_cget_field_id_instance(default_instance_id, comp_name,
      grid_name, field_name);
}

/* ---------------------------------------------------------------------- */

int yac_cget_nbr_comps_instance ( int yac_instance_id ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  return yac_couple_config_get_num_components(couple_config);
}

int yac_cget_nbr_comps ( void ) {
  check_default_instance_id("yac_cget_nbr_comps");
  return yac_cget_nbr_comps_instance(default_instance_id);
}

int yac_cget_nbr_grids_instance ( int yac_instance_id ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  return yac_couple_config_get_num_grids(couple_config);
}

int yac_cget_nbr_grids ( ){
  check_default_instance_id("yac_cget_nbr_grids");
  return yac_cget_nbr_grids_instance( default_instance_id );
}

int yac_cget_comp_nbr_grids_instance ( int yac_instance_id,
                                       const char* comp_name ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  YAC_ASSERT(comp_name != NULL,
    "ERROR(yac_cget_comp_grid_names_instance):"
    "Invalid comp_name. (NULL is not allowed)");
  int nbr_grids = 0;
  for(int i=0; i<yac_couple_config_get_num_grids(couple_config); ++i){
    const char* grid_name = yac_couple_config_get_grid_name(couple_config, i);
    if (yac_cget_nbr_fields_instance(yac_instance_id, comp_name, grid_name) > 0)
      nbr_grids++;
  }
  return nbr_grids;
}

int yac_cget_comp_nbr_grids ( const char* comp_name ){
  check_default_instance_id("yac_cget_comp_nbr_grids");
  return yac_cget_comp_nbr_grids_instance( default_instance_id, comp_name );
}

int yac_cget_nbr_fields_instance ( int yac_instance_id, const char* comp_name,
  const char* grid_name) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);
  int nbr_fields = 0;
  size_t nbr_comp_fields =
    yac_couple_config_get_num_fields(couple_config, comp_idx);
  for(size_t field_idx=0; field_idx<nbr_comp_fields; ++field_idx)
    if(yac_couple_config_field_is_valid(couple_config, comp_idx, field_idx) &&
       !strcmp(
          grid_name,
          yac_couple_config_get_field_grid_name(
            couple_config, comp_idx, field_idx)))
      nbr_fields++;
  return nbr_fields;
}

int yac_cget_nbr_fields ( const char* comp_name, const char* grid_name ) {
  check_default_instance_id("yac_cget_nbr_fields");
  return yac_cget_nbr_fields_instance(default_instance_id, comp_name, grid_name);
}

void yac_cget_comp_names_instance ( int yac_instance_id, int nbr_comps,
  const char ** comp_names ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  for(size_t i=0;i<nbr_comps;++i) {
    comp_names[i] = yac_couple_config_get_component_name(couple_config, i);
  }
}

void yac_cget_comp_names ( int nbr_comps, const char ** comp_names ) {
  check_default_instance_id("yac_cget_comp_names");
  yac_cget_comp_names_instance (default_instance_id, nbr_comps,
    comp_names );
}

void yac_cget_grid_names_instance ( int yac_instance_id,
                                    int nbr_grids, const char ** grid_names ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  for(size_t i=0; i<nbr_grids; ++i){
    grid_names[i] = yac_couple_config_get_grid_name(couple_config, i);
  }
}

void yac_cget_grid_names ( int nbr_grids, const char ** grid_names ) {
  check_default_instance_id("yac_cget_grid_names");
  yac_cget_grid_names_instance(default_instance_id, nbr_grids, grid_names);
}

void yac_cget_comp_grid_names_instance ( int yac_instance_id, const char* comp_name,
                                         int nbr_grids, const char ** grid_names ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  size_t idx = 0;
  YAC_ASSERT(comp_name != NULL,
    "ERROR(yac_cget_comp_grid_names_instance):"
    "Invalid comp_name. (NULL is not allowed)");
  for(size_t i=0; i<yac_couple_config_get_num_grids(couple_config); ++i){
    const char* grid_name = yac_couple_config_get_grid_name(couple_config, i);
    if (yac_cget_nbr_fields_instance(yac_instance_id, comp_name, grid_name) > 0){
      grid_names[idx] = grid_name;
      idx++;
    }
  }
}

void yac_cget_comp_grid_names ( const char* comp_name, int nbr_grids,
                                const char ** grid_names ) {
  check_default_instance_id("yac_cget_grid_names");
  yac_cget_comp_grid_names_instance(default_instance_id, comp_name, nbr_grids, grid_names);
}

void yac_cget_field_names_instance ( int yac_instance_id,
  const char * comp_name, const char* grid_name,
  int nbr_fields, const char ** field_names ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  size_t comp_idx = yac_couple_config_get_component_idx(couple_config, comp_name);

  size_t nbr_comp_fields =
    yac_couple_config_get_num_fields(couple_config, comp_idx);

  for(size_t field_idx=0; field_idx<nbr_comp_fields; ++field_idx) {
    if(yac_couple_config_field_is_valid(couple_config, comp_idx, field_idx) &&
       !strcmp(
          grid_name,
          yac_couple_config_get_field_grid_name(
            couple_config, comp_idx, field_idx))) {
      *field_names =
        yac_couple_config_get_field_name(
          couple_config, comp_idx, field_idx);
      field_names++;
    }
  }
}

void yac_cget_field_names ( const char* comp_name, const char* grid_name,
  int nbr_fields, const char ** field_names ) {
  check_default_instance_id("yac_cget_field_names");
  yac_cget_field_names_instance ( default_instance_id,
    comp_name, grid_name, nbr_fields, field_names );
}

const char* yac_cget_component_name_from_field_id ( int field_id ) {
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  struct component * component = yac_get_coupling_field_component(field);
  return yac_component_get_name(component);
}

const char* yac_cget_grid_name_from_field_id ( int field_id ) {
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  struct yac_basic_grid * grid = yac_coupling_field_get_basic_grid(field);
  return yac_basic_grid_get_name(grid);
}

const char* yac_cget_field_name_from_field_id ( int field_id ) {
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  return yac_get_coupling_field_name(field);
}

const char* yac_cget_timestep_from_field_id ( int field_id ){
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  return yac_get_coupling_field_timestep(field);
}

int yac_cget_collection_size_from_field_id ( int field_id ) {
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  return yac_get_coupling_field_collection_size(field);
}

int yac_cget_role_from_field_id ( int field_id ){
  struct coupling_field * field =
    yac_unique_id_to_pointer(field_id, "field_id");
  YAC_ASSERT(field != NULL, "ERROR: field ID not defined!");
  return yac_get_coupling_field_exchange_type(field);
}

/* ---------------------------------------------------------------------- */

const char* yac_cget_field_timestep_instance(
  int yac_instance_id, const char* comp_name, const char* grid_name,
  const char* field_name ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  return
    yac_couple_config_get_field_timestep(
      couple_config, comp_name, grid_name, field_name);
}

const char* yac_cget_field_timestep ( const char* comp_name, const char* grid_name,
  const char* field_name ) {
  check_default_instance_id("yac_cget_field_timestep");
  return yac_cget_field_timestep_instance(default_instance_id, comp_name,
      grid_name, field_name);
}

double yac_cget_field_frac_mask_fallback_value_instance(
  int yac_instance_id, const char* comp_name, const char* grid_name,
  const char* field_name) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  return
    yac_couple_config_get_frac_mask_fallback_value(
      yac_instance_get_couple_config(instance),
      comp_name, grid_name, field_name);
}

int yac_cget_field_collection_size_instance ( int yac_instance_id,
  const char* comp_name, const char* grid_name, const char* field_name ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  return
    yac_couple_config_get_collection_size(
      yac_instance_get_couple_config(instance),
      comp_name, grid_name, field_name);
}

double yac_cget_field_frac_mask_fallback_value(
  const char* comp_name, const char* grid_name, const char* field_name) {
  check_default_instance_id("yac_cget_field_frac_mask_fallback_value");
  return
    yac_cget_field_frac_mask_fallback_value_instance(
      default_instance_id, comp_name, grid_name, field_name);
}

int yac_cget_field_collection_size ( const char* comp_name,
  const char* grid_name, const char* field_name ) {
  check_default_instance_id("yac_cget_field_collection_size");
  return yac_cget_field_collection_size_instance(default_instance_id, comp_name,
    grid_name, field_name);
}

int yac_cget_field_role_instance ( int yac_instance_id, const char* comp_name,
  const char* grid_name, const char* field_name ) {
  struct yac_instance * instance =
    yac_unique_id_to_pointer(yac_instance_id, "yac_instance_id");
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  return
    yac_couple_config_get_field_role(
      couple_config, comp_name, grid_name, field_name);
}

int yac_cget_field_role ( const char* comp_name, const char* grid_name,
  const char* field_name ) {
  check_default_instance_id("yac_cget_field_role");
  return yac_cget_field_role_instance(default_instance_id, comp_name,
    grid_name, field_name);
}

/* ---------------------------------------------------------------------- */

void yac_cdef_grid_reg2d(
  const char * grid_name, int nbr_vertices[2], int cyclic[2],
  double *x_vertices, double *y_vertices, int *grid_id) {

  size_t nbr_vertices_size_t[2] =
    {(size_t)nbr_vertices[0], (size_t)nbr_vertices[1]};

  check_x_vertices(x_vertices, nbr_vertices[0], "yac_cdef_grid_reg2d");
  check_y_vertices(y_vertices, nbr_vertices[1], "yac_cdef_grid_reg2d");

  *grid_id =
    yac_add_grid(
      grid_name, yac_generate_basic_grid_data_reg_2d(
                   nbr_vertices_size_t, cyclic, x_vertices, y_vertices));
}

/* ---------------------------------------------------------------------- */

void yac_cdef_grid_curve2d(
  const char * grid_name, int nbr_vertices[2], int cyclic[2],
  double *x_vertices, double *y_vertices, int *grid_id) {

  size_t nbr_vertices_size_t[2] =
    {(size_t)nbr_vertices[0], (size_t)nbr_vertices[1]};

  check_x_vertices(x_vertices, nbr_vertices[0], "yac_cdef_grid_curve2d");
  check_y_vertices(y_vertices, nbr_vertices[1], "yac_cdef_grid_curve2d");

  *grid_id =
    yac_add_grid(
      grid_name, yac_generate_basic_grid_data_curve_2d(
                   nbr_vertices_size_t, cyclic, x_vertices, y_vertices));
}

/* ---------------------------------------------------------------------- */

void yac_cdef_grid_unstruct(
  const char * grid_name, int nbr_vertices,
  int nbr_cells, int *num_vertices_per_cell, double *x_vertices,
  double *y_vertices, int *cell_to_vertex, int *grid_id) {

  check_x_vertices(x_vertices, nbr_vertices, "yac_cdef_grid_unstruct");
  check_y_vertices(y_vertices, nbr_vertices, "yac_cdef_grid_unstruct");

  *grid_id =
    yac_add_grid(
      grid_name,
      yac_generate_basic_grid_data_unstruct(
        (size_t)nbr_vertices, (size_t)nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

/* ---------------------------------------------------------------------- */

void yac_cdef_grid_unstruct_ll(
  const char * grid_name, int nbr_vertices,
  int nbr_cells, int *num_vertices_per_cell, double *x_vertices,
  double *y_vertices, int *cell_to_vertex, int *grid_id) {

  check_x_vertices(x_vertices, nbr_vertices, "yac_cdef_grid_unstruct_ll");
  check_y_vertices(y_vertices, nbr_vertices, "yac_cdef_grid_unstruct_ll");

  *grid_id =
    yac_add_grid(
      grid_name,
      yac_generate_basic_grid_data_unstruct_ll(
        (size_t)nbr_vertices, (size_t)nbr_cells, num_vertices_per_cell,
        x_vertices, y_vertices, cell_to_vertex));
}

/* ---------------------------------------------------------------------- */

void yac_cset_global_index(
  int const * global_index, int location, int grid_id) {

  struct basic_grid_data * grid_data =
    yac_basic_grid_get_data(yac_unique_id_to_pointer(grid_id, "grid_id"));

  yac_int ** grid_global_ids;
  size_t count;

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_cset_global_index): invalid location")
  switch (location) {
    case (CELL): {
      grid_global_ids = &(grid_data->cell_ids);
      count = grid_data->num_cells;
      break;
    }
    case (CORNER): {
      grid_global_ids = &(grid_data->vertex_ids);
      count = grid_data->num_vertices;
      break;
    }
    default:
    case (EDGE): {
      grid_global_ids = &(grid_data->edge_ids);
      count = grid_data->num_edges;
      break;
    }
  }

  size_t global_ids_size = count * sizeof(**grid_global_ids);
  if (*grid_global_ids == NULL) *grid_global_ids = xmalloc(global_ids_size);
  for (size_t i = 0; i < count; ++i)
    (*grid_global_ids)[i] = (yac_int)(global_index[i]);
}

/* ---------------------------------------------------------------------- */

void yac_cset_core_mask(
  int const * is_core, int location, int grid_id) {

  struct basic_grid_data * grid_data =
    yac_basic_grid_get_data(yac_unique_id_to_pointer(grid_id, "grid_id"));

  int ** grid_core_mask;
  size_t count;

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_cset_core_mask): invalid location")
  switch (location) {
    case (CELL): {
      grid_core_mask = &(grid_data->core_cell_mask);
      count = grid_data->num_cells;
      break;
    }
    case (CORNER): {
      grid_core_mask = &(grid_data->core_vertex_mask);
      count = grid_data->num_vertices;
      break;
    }
    default:
    case (EDGE): {
      grid_core_mask = &(grid_data->core_edge_mask);
      count = grid_data->num_edges;
      break;
    }
  }

  size_t core_mask_size = count * sizeof(**grid_core_mask);
  if (*grid_core_mask == NULL) *grid_core_mask = xmalloc(core_mask_size);
  memcpy(*grid_core_mask, is_core, core_mask_size);
}

/* ---------------------------------------------------------------------- */

size_t yac_get_grid_size ( int located, int grid_id ) {

  struct yac_basic_grid * grid =
    yac_unique_id_to_pointer(grid_id, "grid_id");
  enum yac_location location = yac_get_location(located);

  return yac_basic_grid_get_data_size(grid, location);
}

/* ---------------------------------------------------------------------- */

size_t yac_get_points_size ( int points_id  ) {

  struct user_input_data_points * points =
    yac_unique_id_to_pointer(points_id, "points_id");
  struct yac_basic_grid * grid = points->grid;
  enum yac_location location = points->location;

  return yac_basic_grid_get_data_size(grid, location);
}

/* ---------------------------------------------------------------------- */

void yac_cget_interp_stack_config(int * interp_stack_config_id) {

  interp_stack_configs =
    xrealloc(
      interp_stack_configs,
      (num_interp_stack_configs + 1) * sizeof(*interp_stack_configs));

  interp_stack_configs[num_interp_stack_configs] =
    yac_interp_stack_config_new();
  *interp_stack_config_id =
    yac_pointer_to_unique_id(interp_stack_configs[num_interp_stack_configs]);
  num_interp_stack_configs++;
}

void yac_cfree_interp_stack_config(int interp_stack_config_id) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");
  yac_interp_stack_config_delete(interp_stack_config);

  int flag = 0;
  for (size_t i = 0; i < num_interp_stack_configs; ++i) {
    if (interp_stack_configs[i] == interp_stack_config)
      interp_stack_configs[i] = NULL;
    flag |= (interp_stack_configs[i] != NULL);
  }
  if (!flag) {
    free(interp_stack_configs);
    interp_stack_configs = NULL;
    num_interp_stack_configs = 0;
  }
}

void yac_cadd_interp_stack_config_average(
  int interp_stack_config_id,
  int reduction_type, int partial_coverage) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  YAC_ASSERT(
    (reduction_type == AVG_ARITHMETIC) ||
    (reduction_type == AVG_DIST) ||
    (reduction_type == AVG_BARY),
    "ERROR(yac_add_interp_stack_config_average): invalid reduction type")

  yac_interp_stack_config_add_average(
    interp_stack_config, (enum yac_interp_avg_weight_type)reduction_type,
    partial_coverage);
}

void yac_cadd_interp_stack_config_nnn(
  int interp_stack_config_id,
  int type, size_t n, double scale) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  YAC_ASSERT(
    (type == NNN_AVG) ||
    (type == NNN_DIST) ||
    (type == NNN_GAUSS) ||
    (type == NNN_RBF),
    "ERROR(yac_add_interp_stack_config_nnn): invalid weightening type")

  yac_interp_stack_config_add_nnn(
    interp_stack_config, (enum yac_interp_nnn_weight_type)type, n, scale);
}

void yac_cadd_interp_stack_config_conservative(
  int interp_stack_config_id, int order, int enforced_conserv,
  int partial_coverage, int normalisation) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  YAC_ASSERT(
    (normalisation == CONSERV_DESTAREA) ||
    (normalisation == CONSERV_FRACAREA),
    "ERROR(yac_add_interp_stack_config_conservative):"
    "invalid normalisation type")

  yac_interp_stack_config_add_conservative(
    interp_stack_config, order, enforced_conserv, partial_coverage,
    (enum yac_interp_method_conserv_normalisation)normalisation);
}

void yac_cadd_interp_stack_config_spmap(
  int interp_stack_config_id, double spread_distance,
  double max_search_distance, int weight_type) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  YAC_ASSERT(
    (weight_type == SPMAP_AVG) ||
    (weight_type == SPMAP_DIST),
    "ERROR(yac_add_interp_stack_config_spmap):"
    "invalid weightening type")

  yac_interp_stack_config_add_spmap(
    interp_stack_config, spread_distance, max_search_distance,
    (enum yac_interp_spmap_weight_type)weight_type);
}

void yac_cadd_interp_stack_config_hcsbb(int interp_stack_config_id) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_hcsbb(interp_stack_config);
}

void yac_cadd_interp_stack_config_user_file(
  int interp_stack_config_id, char const * filename, char const * src_grid_name,
  char const * tgt_grid_name) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_user_file(
    interp_stack_config, filename, src_grid_name,
    tgt_grid_name);
}

void yac_cadd_interp_stack_config_fixed(
  int interp_stack_config_id, double value) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_fixed(interp_stack_config, value);
}

void yac_cadd_interp_stack_config_check(
  int interp_stack_config_id, char const * constructor_key,
  char const * do_search_key) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_check(
    interp_stack_config, constructor_key, do_search_key);
}

void yac_cadd_interp_stack_config_creep(
  int interp_stack_config_id, int creep_distance) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_creep(interp_stack_config, creep_distance);
}

void yac_cadd_interp_stack_config_user_callback(
  int interp_stack_config_id, char const * func_compute_weights_key) {

  struct yac_interp_stack_config * interp_stack_config =
    yac_unique_id_to_pointer(
      interp_stack_config_id, "interp_stack_config_id");

  yac_interp_stack_config_add_user_callback(
    interp_stack_config, func_compute_weights_key);
}

/* ---------------------------------------------------------------------- */

void yac_cadd_compute_weights_callback(
  yac_func_compute_weights compute_weights_callback,
  void * user_data, char const * key) {

  yac_interp_method_callback_add_compute_weights_callback(
    compute_weights_callback, user_data, key);
}
