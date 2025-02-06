// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONFIG_YAML_H
#define CONFIG_YAML_H

#include "couple_config.h"

extern int const YAC_YAML_PARSER_DEFAULT;    //!<  default parse flags (YAML format)
extern int const YAC_YAML_PARSER_JSON_AUTO;  //!<  switch to JSON format,
                                             //!<  if indicated by file extension
extern int const YAC_YAML_PARSER_JSON_FORCE; //!<  assume JSON format

/**
 * Reader for yaml while, which parses the given yaml file and adds the
 * coupling from the file to the coupling configuration data
 *
 * @param[in,out] couple_config coupling configuration data
 * @param[in]     yaml_filename name of yaml configuration file
 * @param[in]     parse_flags   flags to be used for parsing the
 *                              configuration file
 */
void yac_yaml_read_coupling(
  struct yac_couple_config * couple_config, const char * yaml_filename,
  int parse_flags);

/**
 * Emit coupling configuration to string
 *
 * @param[in] couple_config coupling configuration
 * @param[in] emit_flags    flags to be used for emitting the
 *                          coupling configuration
 * @param[in] include_definitions include user definitions (components, grids,
 *                                and fields) in the output file
 * @return string containing coupling configuration
 */
char * yac_yaml_emit_coupling(
  struct yac_couple_config * couple_config, int emit_flags,
  int include_definitions);

/**
 * Parse a "0"-terminated string that contains a interpolation stack
 * configuration
 *
 * @param[in] interp_stack_config string containing interpolation stack
 *                                configuration
 * @param[in] parse_flags         flags to be used for parsing the
 *                                string
 */
struct yac_interp_stack_config *
  yac_yaml_parse_interp_stack_config_string(
    char const * interp_stack_config, int parse_flags);

#endif /* CONFIG_YAML_H */

