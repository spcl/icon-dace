// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// YAC PUBLIC HEADER START

/**
 * reads a coordinate array from a netcdf file, tries to interpret the units
 * attribute of the variable and will convert the coordinates into radiant
 * if necessary
 * @param[in]  ncid
 * @param[in]  lon_name name of longitude array in grid file
 * @param[in]  lat_name name of latitude array in grid file
 * @param[out] lon      lon array in radiant
 * @param[out] lat      lat array in radiant
 * @param[out] len      number of coordinates in lon and lat arrays
 */
void yac_read_coords(
  int ncid, char const * lon_name, char const * lat_name,
  double ** lon, double ** lat, size_t * len);

/**
 * Checks the variable for an attribute with the name "units". If it is
 * available and contains the string "degree", this routine will return 1,
 * 0 otherwise.
 */
int yac_check_coord_units(int ncid, int varid);

// YAC PUBLIC HEADER STOP
