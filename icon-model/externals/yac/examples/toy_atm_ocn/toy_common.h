// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

const char * fieldName[] = {"conserv_field_atm_oce",
                            "hcsbb_field_atm_oce",
                            "conserv_field_oce_atm",
                            "hcsbb_field_oce_atm"};
int const num_fields = (int)(sizeof(fieldName)/sizeof(fieldName[0]));
const int num_steps = 1;

static inline void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

static inline void XYZtoLL (double const p_in[], double * lon, double * lat) {

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))
