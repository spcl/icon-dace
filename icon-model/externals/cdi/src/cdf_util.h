#ifndef CDF_UTIL_H_
#define CDF_UTIL_H_

#include <stdbool.h>

bool xtypeIsText(int xtype);

int get_time_units(size_t len, const char *ptu);

bool is_time_units(const char *timeunits);
bool is_timeaxis_units(const char *timeunits);

bool is_height_units(const char *units);
bool is_pressure_units(const char *units);
bool is_DBL_axis(/*const char *units,*/ const char *longname);
bool is_depth_axis(const char *stdname, const char *longname);
bool is_height_axis(const char *stdname, const char *longname);
bool is_altitude_axis(const char *stdname, const char *longname);
bool is_reference_axis(const char *stdname, const char *longname);

bool is_lon_axis(const char *units, const char *stdname);
bool is_lat_axis(const char *units, const char *stdname);

bool is_x_axis(const char *units, const char *stdname);
bool is_y_axis(const char *units, const char *stdname);

void cdf_set_gridtype(const char *attstring, int *gridtype);
void cdf_set_zaxistype(const char *attstring, int *zaxistype);
int attribute_to_calendar(const char *attstring);

#endif
