"""
This module provides a function to copy a NetCDF file while converting variable
data types.
"""

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from util.log_handler import logger


def nc4_get_copy(name_in, name_out):
    din = Dataset(name_in, "r")
    logger.info("creating netCDF4 file: %s", name_out)
    d_out = Dataset(name_out, "w")

    # Copy dimensions
    for dname, dim in din.dimensions.items():
        d_out.createDimension(dname, len(dim) if not dim.isunlimited() else None)

    # Copy variables
    for v_name, in_var in din.variables.items():
        # always create double precision output.
        # Default is single which will not work for perturb.
        out_var = d_out.createVariable(v_name, np.float64, in_var.dimensions)

        # Copy variable attributes, make sure float32 attributes are created as float64
        real_atts = [
            k for k in in_var.ncattrs() if isinstance(in_var.getncattr(k), np.float32)
        ]
        ncdict = {
            k: in_var.getncattr(k) for k in in_var.ncattrs() if k not in real_atts
        }
        ncdict_dbl = {
            k: np.float64(in_var.getncattr(k))
            for k in in_var.ncattrs()
            if k in real_atts
        }
        out_var.setncatts({**ncdict, **ncdict_dbl})

        out_var[:] = np.float64(in_var[:])

    # close the input file
    din.close()
    return d_out
