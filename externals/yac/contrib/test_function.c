/**
 * @file test_function.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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

#ifndef TEST_FUNCTION
#define TEST_FUNCTION
#include <math.h>
#include "grid.h"
#include "geometry.h"

double test_func(double lon, double lat) {

  double tmp = sin(lon * 2.0);
  return tmp * cos(lat * 3.0);
}

double test_func_deg(double lon, double lat) {

  return test_func(lon * YAC_RAD, lat * YAC_RAD);
}

double test_ana_fcos(double lon, double lat) {

  double dp_length = 1.2 * M_PI;
  double coef = 2.0;
  double coefmult = 1.0;
  return coefmult * (coef - cos(M_PI * acos(cos(lon)*cos(lat) )/dp_length));
}

double test_ana_fcossin(double lon, double lat) {

  double dp_length = 1.0 * M_PI;
  double coef = 21.0;
  double coefmult = 3.846 * 20.0;
  return coefmult * (coef - cos(M_PI*(acos(cos(lon)*cos(lat))/dp_length)) *
                            sin(M_PI*(asin(sin(lon)*sin(lat))/dp_length)));
}

double test_one(double lon, double lat) {
  return 1.0;
}

double test_gulfstream(double lon, double lat) {

  // Analytical Gulf Stream
  double gf_coef = 1.0; // Coefficient for Gulf Stream term (0.0 = no Gulf Stream)
  double gf_ori_lon = -80.0 * YAC_RAD; // Origin of the Gulf Stream
  double gf_ori_lat =  25.0 * YAC_RAD; // Origin of the Gulf Stream
  double gf_end_lon =  -1.8 * YAC_RAD; // End point of the Gulf Stream
  double gf_end_lat =  50.0 * YAC_RAD; // End point of the Gulf Stream
  double gf_dmp_lon = -25.5 * YAC_RAD; // Point of the Gulf Stream decrease
  double gf_dmp_lat =  55.5 * YAC_RAD; // Point of the Gulf Stream decrease

  double dr0 = sqrt(pow(gf_end_lon - gf_ori_lon, 2.0) +
                    pow(gf_end_lat - gf_ori_lat, 2.0));
  double dr1 = sqrt(pow(gf_dmp_lon - gf_ori_lon, 2.0) +
                    pow(gf_dmp_lat - gf_ori_lat, 2.0));

  if (lon > M_PI) lon -= 2.0 * M_PI;
  if (lon < -M_PI) lon += 2.0 * M_PI;

  double dx = lon - gf_ori_lon;
  double dy = lat - gf_ori_lat;
  double dr = sqrt(dx * dx + dy * dy);
  double dth = atan2(dy, dx);
  double dc = 1.3 * gf_coef;
  if (dr > dr0) dc = 0.0;
  if (dr > dr1) dc = dc * cos(M_PI*0.5*(dr-dr1)/(dr0-dr1));
  return test_ana_fcos(lon, lat) +
         (MAX(1000.0 * sin(0.4 * (0.5*dr+dth) +
          0.007*cos(50.0*dth) + 0.37*M_PI),999.0) - 999.0) * dc;
}

double test_harmonic(double lon, double lat) {
  return 2.0 + pow(sin(2.0 * lat), 16.0)  * cos(16.0 * lon);
}

double test_vortex(double lon, double lat) {
  double dLon0 = 5.5;
  double dLat0 = 0.2;
  double dR0   = 3.0;
  double dD    = 5.0;
  double dT    = 6.0;

  double dSinC = sin(dLat0);
  double dCosC = cos(dLat0);

  // Find the rotated longitude and latitude of a point on a sphere
  //		with pole at (dLon0, dLat0).
  double dCosT = cos(lat);
  double dSinT = sin(lat);

  double dTrm = dCosT * cos(lon - dLon0);
  double dX   = dSinC * dTrm - dCosC * dSinT;
  double dY   = dCosT * sin(lon - dLon0);
  double dZ   = dSinC * dSinT + dCosC * dTrm;

  double dlon = atan2(dY, dX);
  if(dlon < 0.0) dlon = dlon + 2.0 * M_PI;
  double dlat = asin(dZ);

  double dRho = dR0 * cos(dlat);
  double dVt = 3.0 * sqrt(3.0)/2.0/cosh(dRho)/cosh(dRho)*tanh(dRho);
  double dOmega = (dRho == 0.0)?0.0:dVt/dRho;

  return 2.0 * (1.0 + tanh(dRho / dD * sin(dlon - dOmega * dT)));
}

#endif // TEST_FUNCTION

