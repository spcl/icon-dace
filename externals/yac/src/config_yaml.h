/**
 * @file config_yaml.h
 *
 * Interface for retrieving values given in coupling configuration YAML file.
 *
 * @copyright Copyright  (C)  2022 Moritz Hanke <hanke@dkrz.de>
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

#ifndef CONFIG_YAML_H
#define CONFIG_YAML_H


#include "couple_config.h"

extern int YAC_YAML_PARSER_DEFAULT;    //!<  default parse flags (YAML format)
extern int YAC_YAML_PARSER_JSON_AUTO;  //!<  switch to JSON format,
                                       //!<  if indicated by file extension
extern int YAC_YAML_PARSER_JSON_FORCE; //!<  assume JSON format

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
 * @return string containing coupling configuration
 */
char * yac_yaml_emit_coupling(
  struct yac_couple_config * couple_config, int emit_flags);

#endif /* CONFIG_YAML_H */

