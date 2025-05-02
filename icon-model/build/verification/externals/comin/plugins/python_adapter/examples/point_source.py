"""
Test plugin for the ICON Community Interface (ComIn)

The point_source plugin shows an example for
- Requesting a tracer that participates in ICON's turbulence and convection scheme
- Adding point source emissions to this tracer
- Using KDTree of scipy to locate the point source
- Use decomp_domain to check if the point source is local to the PE
- Updating the tracer with tendencies received from ICON's turbulence and convection scheme

@authors 12/2023 :: ICON Community Interface  <comin@icon-model.org>

SPDX-License-Identifier: BSD-3-Clause

Please see the file LICENSE in the root of the source tree for this code.
Where software is supplied by third parties, it is indicated in the
headers of the routines.
"""

import comin
import numpy as np
import sys
from scipy.spatial import KDTree

# request to add point source tracer on all domains
comin.var_request_add(("pntsrc", -1), False)
comin.metadata_set(("pntsrc", -1), tracer=True, tracer_turb=True, tracer_conv=True)

# Point source is only added to domain 1
jg = 1
msgrank = 0 # Rank that prints messages
pntsrc_lon = 141.0325
pntsrc_lat = 37.421389

def message(message_string, rank):
    """Short helper function to print a message on one PE"""
    if (comin.parallel_get_host_mpi_rank() == rank):
        print(f"ComIn point_source.py: {message_string}", file=sys.stderr)

def lonlat2xyz(lon, lat):
    clat = np.cos(lat)
    return clat * np.cos(lon), clat * np.sin(lon), np.sin(lat)

@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def pntsrc_constructor():
    """Constructor: Get pointers to tracers and descriptive data"""
    global pntsrc, ddt_pntsrc_turb, ddt_pntsrc_conv

    entry_points = [comin.EP_ATM_ADVECTION_BEFORE,comin.EP_ATM_PHYSICS_BEFORE]
    pntsrc = comin.var_get(entry_points, ("pntsrc", jg),
                           comin.COMIN_FLAG_READ | comin.COMIN_FLAG_WRITE)
    # ICON prepends 'ddt_' and appends '_turb'/'_conv' for the respective tendencies
    entry_points = [comin.EP_ATM_PHYSICS_BEFORE]
    ddt_pntsrc_turb = comin.var_get(entry_points,("ddt_pntsrc_turb", jg), comin.COMIN_FLAG_READ)
    ddt_pntsrc_conv = comin.var_get(entry_points,("ddt_pntsrc_conv", jg), comin.COMIN_FLAG_READ)

    message("pntsrc_constructor successful", msgrank)

@comin.register_callback(comin.EP_ATM_INIT_FINALIZE)
def pntsrc_init():
    global lcontain_pntsrc, jc_loc, jb_loc

    # pntsrc_init is executed outside domain loop
    # all arrays are for domain 1 only
    domain = comin.descrdata_get_domain(jg)
    clon = np.asarray(domain.cells.clon)
    clat = np.asarray(domain.cells.clat)
    xyz = np.c_[lonlat2xyz(clon.ravel(),clat.ravel())]
    decomp_domain = np.asarray(domain.cells.decomp_domain)

    tree = KDTree(xyz)
    dd, ii = tree.query([lonlat2xyz(np.deg2rad(pntsrc_lon), np.deg2rad(pntsrc_lat))], k=1)

    lcontain_pntsrc = False
    if (decomp_domain.ravel()[ii] == 0):
        # point found is inside prognostic area
        # This implicitly assumes that on each other PE, the nearest neighbor is located in the halo zone
        jc_loc, jb_loc = np.unravel_index(ii, clon.shape)
        message(f"Min at PE {comin.parallel_get_host_mpi_rank()}, clon={np.rad2deg(clon[jc_loc,jb_loc])}, clat={np.rad2deg(clat[jc_loc,jb_loc])}", comin.parallel_get_host_mpi_rank())
        lcontain_pntsrc = True

@comin.register_callback(comin.EP_ATM_ADVECTION_BEFORE)
def pntsrc_emission():
    """Emission for pntsrc on domain 1"""
    if (comin.current_get_domain_id() == jg and lcontain_pntsrc):
        pntsrc_np = np.asarray(pntsrc)
        pntsrc_np[jc_loc,:,jb_loc,0,0] = pntsrc_np[jc_loc,:,jb_loc,0,0] + 1.

@comin.register_callback(comin.EP_ATM_PHYSICS_BEFORE)
def pntsrc_updatePhysTends():
    """Update tracer from turb and conv tendencies for domain 1"""
    if (comin.current_get_domain_id() == jg):
        pntsrc_np = np.asarray(pntsrc)
        ddt_pntsrc_turb_np = np.asarray(ddt_pntsrc_turb)
        ddt_pntsrc_conv_np = np.asarray(ddt_pntsrc_conv)

        dtime = comin.descrdata_get_timesteplength(jg)
        pntsrc_np[:,:,:,0,0] = np.maximum( (pntsrc_np[:,:,:,0,0] \
            + dtime * (ddt_pntsrc_turb_np[:,:,:,0,0] + ddt_pntsrc_conv_np[:,:,:,0,0])), 0.)

@comin.register_callback(comin.EP_DESTRUCTOR)
def pntsrc_destructor():
    message("pntsrc_destructor called!", msgrank)
