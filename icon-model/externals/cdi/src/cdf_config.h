#ifndef CDF_CONFIG_H_
#define CDF_CONFIG_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include <netcdf.h>

#ifdef NC_FORMAT_64BIT_DATA
#define HAVE_NETCDF5 1
#endif

#endif

#endif
