/**
 * @file test_cxc.h
 *
 * @copyright Copyright  (C)  2014 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifndef TEST_CXC_H
#define TEST_CXC_H

#include "geometry.h"
#include "tests.h"

/** \example test_cxc.c
 * A test for cell intersections.
 */

void test_cxc(double lon_a, double lat_a, double lon_b, double lat_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
              double lon_ref_p, double lat_ref_p,
              double lon_ref_q, double lat_ref_q, int ref_ret_val);

void test_cxc_rad(double lon_a, double lat_a, double lon_b, double lat_b,
                  double lon_c, double lat_c, double lon_d, double lat_d,
                  enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
                  double lon_ref_p, double lat_ref_p,
                  double lon_ref_q, double lat_ref_q, int ref_ret_val);
void get_edge_middle_point(
  enum yac_edge_type edge_type,
  double lon_a, double lat_a, double lon_b, double lat_b,
  double * lon_middle, double * lat_middle);

#endif // TEST_CXC_H

