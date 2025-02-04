#!/usr/bin/env python
import os
from typing import Dict, List, Tuple
from collections import defaultdict
from pathlib import Path
import argparse

from yaml import dump as dump_yaml, load as load_yaml
try:
    from yaml import CDumper as YAML_Dumper, CLoader as YAML_Loader
except ImportError:
    from yaml import Dumper as YAML_Dumper, Loader as YAML_Loader

import dace


###############################################################################
# Start demo SDFG (edges_to_cells_bilinear_interpolation)
###############################################################################
# Original ICON code given below. There are two major differences in the DaCe Python version:
# 1. We have to change Fortran 1-based indexing to DaCe 0-based indexing
# 2. Fortran shapes and indices are reversed in DaCe (column vs row major ordering)
#
#   #ifdef __LOOP_EXCHANGE
#         DO jc = i_startidx, i_endidx
#           DO jk = 1, nlev
#           z_ekinh(jk,jc,jb) =  &
#   #else
#         DO jk = 1, nlev
#           DO jc = i_startidx, i_endidx
#           z_ekinh(jc,jk,jb) =  &
#   #endif
#             p_int%e_bln_c_s(jc,1,jb)*z_kin_hor_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
#             p_int%e_bln_c_s(jc,2,jb)*z_kin_hor_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
#             p_int%e_bln_c_s(jc,3,jb)*z_kin_hor_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))
#           ENDDO
#         ENDDO

MIXED_PRECISION = False # `False` for most relevant configurations
LOOP_EXCHANGE = False # usually `True` for CPU, and `False` for GPU

wp = dace.float64 # working precision = double precision
vp = dace.float32 if MIXED_PRECISION else dace.float64 # variable precision

block_size = dace.symbol("block_size")
# usually `num_cell_blocks == 1` for GPU configurations
num_cell_blocks = dace.symbol("num_cell_blocks")
# usually `num_edge_blocks == ceil(1.5) = 2` for GPU configurations
num_edge_blocks = dace.symbol("num_edge_blocks")
num_levels = dace.symbol("num_levels")

cell_block_idx = dace.symbol("cell_block_idx")

InterpolationData = dace.data.Structure(
    # Name has to be the same as ICON name, so integration works out
    name="t_int_state",
    members=dict(
        # coefficient for bilinear interpolation from edges to cells for scalar quantities
        # Name has to be the same as ICON name, so integration works out
        e_bln_c_s=wp[num_cell_blocks, 3, block_size],
    ),
)

LOOP_EXCHANGED_VP_CELL_FIELD_TYPE = (
    vp[num_cell_blocks, block_size, num_levels]
    if LOOP_EXCHANGE else
    vp[num_cell_blocks, num_levels, block_size]
)


@dace.program
def edges_to_cells_bilinear_interpolation(
    cell_scalar: LOOP_EXCHANGED_VP_CELL_FIELD_TYPE,
    edge_scalar: vp[num_edge_blocks, num_levels, block_size],
    # lookup tables for local neighbors (3 edges for each triangle/cell)
    # (Note: they are 1-based indices, so we change them to 0-based in the code)
    cell_to_edges_idx: dace.int32[3, num_cell_blocks, block_size],
    cell_to_edges_block_idx: dace.int32[3, num_cell_blocks, block_size],
    # Struct with interpolation coefficients
    interpolation_data: InterpolationData,
    # where to compute
    cell_start_idx: dace.int32,
    cell_end_idx: dace.int32,
    # for simplicity, this stencil only works for one single block
    # (though for GPU configurations, there is only one block encompassing all cells)
    cell_block_idx: dace.int32,
):
    # Change from Fortran 1-based indexing to DaCe 0-based indexing
    cell_start_idx -= 1
    cell_end_idx -= 1
    cell_block_idx -= 1

    # In ICON, order of nesting is determined by `LOOP_EXCHANGE`
    for cell_idx, level in dace.map[cell_start_idx:cell_end_idx + 1, 0:num_levels]:

        interpolated_cell_value = (
            # edge 1
            interpolation_data.e_bln_c_s[cell_block_idx, 0, cell_idx] *
            edge_scalar[
                cell_to_edges_block_idx[0, cell_block_idx, cell_idx] - 1,
                level,
                cell_to_edges_idx[0, cell_block_idx, cell_idx] - 1,
            ]
            +
            # edge 2
            interpolation_data.e_bln_c_s[cell_block_idx, 1, cell_idx] *
            edge_scalar[
                cell_to_edges_block_idx[1, cell_block_idx, cell_idx] - 1,
                level,
                cell_to_edges_idx[1, cell_block_idx, cell_idx] - 1,
            ]
            +
            # edge 3
            interpolation_data.e_bln_c_s[cell_block_idx, 2, cell_idx] *
            edge_scalar[
                cell_to_edges_block_idx[2, cell_block_idx, cell_idx] - 1,
                level,
                cell_to_edges_idx[2, cell_block_idx, cell_idx] - 1,
            ]
        )

        if LOOP_EXCHANGE:
            cell_scalar[cell_block_idx, cell_idx, level] = interpolated_cell_value
        else:
            cell_scalar[cell_block_idx, level, cell_idx] = interpolated_cell_value


EDGES_TO_CELLS_BILINEAR_INTERPOLATION_MODULE_DEFINITIONS = dict(
    # helpers for fortan interface
    warning = "mo_exception",
    MAX_CHAR_LENGTH = "mo_impl_constants",
    # The only symbol needed for the actual SDFG
    t_int_state = "mo_intp_data_strc",
)


EDGES_TO_CELLS_BILINEAR_INTERPOLATION_ASSOCIATIONS = {
    "mo_velocity_advection.f90": {
        432: {
            448: dict(
                block_size = "nproma",
                num_cell_blocks = "p_patch%nblks_c",
                num_levels = "nlev",
                cell_scalar = "z_ekinh",
                edge_scalar = "z_kin_hor_e",
                cell_to_edges_idx = "ieidx",
                cell_to_edges_block_idx = "ieblk",
                interpolation_data = "p_int",
                cell_start_idx = "i_startidx",
                cell_end_idx = "i_endidx",
                cell_block_idx = "jb",
            )
        }
    }
}


def create_edges_to_cells_bilinear_interpolation_sdfg(
    output_folder: Path,
    associations: Dict[str, List[Tuple[int, int]]],
):

    if associations != {"src/atm_dyn_iconam/mo_velocity_advection.f90": [(432, 448)]}:
        # hard-coded associations only works for this integration
        raise NotImplementedError("Demo SDFG doesn't support this integration!")

    # create SDFG
    demo_sdfg = edges_to_cells_bilinear_interpolation.to_sdfg()
    demo_sdfg.save(str(output_folder / "edges_to_cells_bilinear_interpolation_unsimplified.sdfgz"))

    # create module definitions
    with open(output_folder / "edges_to_cells_bilinear_interpolation_module_definitions.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            EDGES_TO_CELLS_BILINEAR_INTERPOLATION_MODULE_DEFINITIONS,
            Dumper=YAML_Dumper,
        ))

    # create associations
    with open(output_folder / "edges_to_cells_bilinear_interpolation_associations.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            EDGES_TO_CELLS_BILINEAR_INTERPOLATION_ASSOCIATIONS,
            Dumper=YAML_Dumper,
        ))

def create_get_albedos_sdfg(
    output_folder: Path,
    associations: Dict[str, List[Tuple[int, int]]],
):
    # create SDFG
    current_folder = os.path.dirname(os.path.abspath(__file__))
    demo_sdfg = dace.SDFG.from_file(current_folder + "/get_albedos.sdfgz")
    demo_sdfg.save(str(output_folder / "get_albedos_unsimplified.sdfgz"))

    # create module definitions
    with open(output_folder / "get_albedos_module_definitions.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            GET_ALBEDOS_MODULE_DEFINITIONS,
            Dumper=YAML_Dumper,
        ))

    # create associations
    with open(output_folder / "get_albedos_associations.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            GET_ALBEDOS_ASSOCIATIONS,
            Dumper=YAML_Dumper,
        ))


VELOCITY_TENDENCIES_MODULE_DEFINITIONS = dict(
    # The only symbol needed for the actual SDFG
    t_patch = "mo_model_domain",
    t_grid_cells = "mo_model_domain",
    t_grid_edges = "mo_model_domain",
    t_grid_vertices = "mo_model_domain",
    t_subset_range = "mo_model_domain",
    t_tangent_vectors = "mo_model_domain",

    t_grid_geometry_info = "mo_grid_geometry_info",

    t_geographical_coordinates = "mo_math_types",
    t_cartesian_coordinates = "mo_math_types",

    t_comm_pattern = "mo_communication",
    t_comm_gather_pattern = "mo_communication",
    t_scatterPattern = "mo_communication",
    t_comm_pattern_collection = "mo_communication",
    t_p_comm_pattern = "mo_communication",

    t_uuid = "mo_util_uuid_types",

    t_grid_domain_decomp_info = "mo_decomposition_tools",
    t_glb2loc_index_lookup = "mo_decomposition_tools",

    t_distrib_read_data = "mo_read_netcdf_types",

    t_dist_dir = "mo_dist_dir",

    t_int_state = "mo_intp_data_strc",
    t_lsq = "mo_intp_data_strc",
    t_cell_environ = "mo_intp_data_strc",
    t_gauss_quad = "mo_intp_data_strc",

    t_nh_state = "mo_nonhydro_types",
    t_nh_diag = "mo_nonhydro_types",
    t_nh_metrics = "mo_nonhydro_types",
    t_nh_prog = "mo_nonhydro_types",
    t_nh_ref = "mo_nonhydro_types",

    t_prepare_adv = "mo_prepadv_types",

    # helpers for fortan interface
    message = "mo_exception",
    em_error = "mo_exception",
    MAX_CHAR_LENGTH = "mo_impl_constants",

    # globals
    dp = "mo_kind",
    wp = "mo_kind",
    sp = "mo_kind",
    vp = "mo_kind",
    i4 = "mo_kind",

    em_none = "mo_exception",
    em_info = "mo_exception",
    em_warn = "mo_exception",
    em_param = "mo_exception",
    em_debug = "mo_exception",
    number_of_warnings = "mo_exception",
    number_of_errors = "mo_exception",
    warning = "mo_exception",

    grid_length_rescale_factor = "mo_grid_config",
    grid_rescale_factor = "mo_grid_config",
    grid_angular_velocity = "mo_grid_config",
    grid_sphere_radius = "mo_grid_config",
    lrescale_ang_vel = "mo_grid_config",
    namelist_grid_angular_velocity = "mo_grid_config",
    no_of_dynamics_grids = "mo_grid_config",
    n_dom = "mo_grid_config",
    n_dom_start = "mo_grid_config",

    # they are private
    #t0sl_bg = "mo_vertical_grid",
    #h_scal_bg = "mo_vertical_grid",
    #del_t_bg = "mo_vertical_grid",

    sphere_geometry = "mo_grid_geometry_info",
    triangular_cell = "mo_grid_geometry_info",

    dbl_eps = "mo_math_constants",

    rd = "mo_physical_constants",
    cpd = "mo_physical_constants",
    cvd = "mo_physical_constants",
    cvd_o_rd = "mo_physical_constants",
    earth_radius = "mo_physical_constants",
    earth_angular_velocity = "mo_physical_constants",
    alf = "mo_physical_constants",
    als = "mo_physical_constants",
    alsdcp = "mo_physical_constants",
    alv = "mo_physical_constants",
    alvdcp = "mo_physical_constants",
    clw = "mo_physical_constants",
    con0_h = "mo_physical_constants",
    con_h = "mo_physical_constants",
    con_m = "mo_physical_constants",
    cpv = "mo_physical_constants",
    cv_i = "mo_physical_constants",
    cv_v = "mo_physical_constants",
    cvv = "mo_physical_constants",
    dv0 = "mo_physical_constants",
    eta0d = "mo_physical_constants",
    grav = "mo_physical_constants",
    o_m_rdv = "mo_physical_constants",
    p0ref = "mo_physical_constants",
    rcpd = "mo_physical_constants",
    rcpl = "mo_physical_constants",
    rcpv = "mo_physical_constants",
    rcvd = "mo_physical_constants",
    rd_o_cpd = "mo_physical_constants",
    rdv = "mo_physical_constants",
    rhoh2o = "mo_physical_constants",
    rhoice = "mo_physical_constants",
    rv = "mo_physical_constants",
    t3 = "mo_physical_constants",
    tmelt = "mo_physical_constants",
    vtmpc1 = "mo_physical_constants",
    vtmpc2 = "mo_physical_constants",

    filename_max = "mo_io_units",
    nerr = "mo_io_units",
    nlog = "mo_io_units",

    mpi_comm_null = "mo_mpi",
    p_mpi_comm_null = "mo_mpi",
    comm_lev = "mo_mpi",
    p_io = "mo_mpi",
    p_pe = "mo_mpi",
    proc_split = "mo_mpi",
    mpi_request_null = "mo_mpi",
    process_mpi_stdio_id = "mo_mpi",
    i_am_accel_node = "mo_mpi",

    max_hw = "mo_impl_constants",
    min_rlcell_int = "mo_impl_constants",
    min_rlcell = "mo_impl_constants",
    min_rlvert_int = "mo_impl_constants",
    min_rlvert = "mo_impl_constants",
    min_rledge_int = "mo_impl_constants",
    min_rledge = "mo_impl_constants",
    max_dom = "mo_impl_constants",
    rayleigh_classic = "mo_impl_constants",
    rayleigh_klemp = "mo_impl_constants",

    grf_bdywidth_c = "mo_impl_constants_grf",
    grf_bdywidth_e = "mo_impl_constants_grf",
    grf_bdywidth_v = "mo_impl_constants_grf",

    nproma = "mo_parallel_config",
    nblocks_c = "mo_parallel_config",
    ignore_nproma_use_nblocks_c = "mo_parallel_config",
    p_test_run = "mo_parallel_config",
    l_test_openmp = "mo_parallel_config",
    icon_comm_openmp = "mo_parallel_config",
    num_io_procs = "mo_parallel_config",
    num_restart_procs = "mo_parallel_config",
    num_prefetch_proc = "mo_parallel_config",
    proc0_shift = "mo_parallel_config",
    use_dycore_barrier = "mo_parallel_config",

    is_iau_active = "mo_initicon_config",
    iau_wgt_dyn = "mo_initicon_config",

    edge2cell_coeff_cc = "mo_intp_data_strc",
)

GET_ALBEDOS_MODULE_DEFINITIONS = VELOCITY_TENDENCIES_MODULE_DEFINITIONS

sub_dict = dict(
    p_prog="p_nh%prog(nnew)",
    p_patch="p_patch" ,
    p_int="p_int",
    p_metrics="p_nh%metrics",
    p_diag="p_nh%diag",
    z_w_concorr_me="z_w_concorr_me",
    z_kin_hor_e="z_kin_hor_e",
    z_vt_ie="z_vt_ie",
    ntnd="ntl1",
    istep="istep",
    lvn_only="transfer(lvn_only, mold=int(1, kind=4))",
    dtime="dtime",
    dt_linintp_ubc="dt_linintp_ubc_nnew",
    lextra_diffu="1",
    lvert_nest="transfer(lvert_nest, mold=int(1, kind=4))",

)
sub_dict2 = dict(
    p_prog="p_nh%prog(nnew)",
    p_patch="p_patch" ,
    p_int="p_int",
    p_metrics="p_nh%metrics",
    p_diag="p_nh%diag",
    z_w_concorr_me="z_w_concorr_me",
    z_kin_hor_e="z_kin_hor_e",
    z_vt_ie="z_vt_ie",
    ntnd="ntl2",
    istep="istep",
    lvn_only="transfer(lvn_only, mold=int(1, kind=4))",
    dtime="dtime",
    dt_linintp_ubc="dt_linintp_ubc_nnew",
    lextra_diffu="1",
    lvert_nest="transfer(lvert_nest, mold=int(1, kind=4))",
)
"""
get_albedos(this_var_378, istartcol_var_380, iendcol_var_381, 
            config_var_379, sw_albedo_direct_var_383, 
            sw_albedo_diffuse_var_384, lw_albedo_var_382)

    call single_level%get_albedos(istartcol, iendcol, config, &
    &                        sw_albedo_direct, sw_albedo_diffuse, &
    &                        lw_albedo)
    to:

"""
sub_dict3 = dict(
    this_var_278="single_level",
    istartcol_var_380="istartcol",
    iendcol_var_381="iendcol",
    config_var_379="config",
    sw_albedo_direct_var_383="sw_albedo_direct",
    sw_albedo_diffuse_var_384="sw_albedo_diffuse",
    lw_albedo_var_382="lw_albedo",
)
"""
p_diag="p_nh%diag",
p_int="p_int",
p_metrics="p_nh%metrics",
p_patch="p_patch",
p_prog="p_nh%prog(nnow)",
z_kin_hor_e="z_kin_hor_e",
z_vt_ie="z_vt_ie",
z_w_concorr_me="z_w_concorr_me",
dt_linintp_ubc="dt_linintp_ubc_nnow",
dtime="dtime",
istep="istep",
lvn_only="transfer(lvn_only, mold=int(1, kind=4))",
ntnd="ntl1",
lextra_diffu="1", # checked, value is good
lvert_nest="transfer(lvert_nest, mold=int(1, kind=4))",
# nproma="nproma",
# nproma_0="nproma",
# nproma_1="nproma",
# nproma_2="nproma",
# nproma_3="nproma",
# nproma_4="nproma",
#edge2cell_coeff_cc="ptr_int%edge2cell_coeff_cc"
"""
VELOCITY_TENDENCIES_ASSOCIATIONS = {
    "mo_solve_nonhydro.f90": {
        465: {
            466: sub_dict
        },
        497: {
            498: sub_dict2
        }
    }
}

GET_ALBEDOS_ASSOCIATIONS = {
    "radiation_interface.F90": {
        335: {
            338: sub_dict3
        },
    }
}

def create_velocity_tendencies_sdfg(
    output_folder: Path,
    associations: Dict[str, List[Tuple[int, int]]],
):

    if associations != {"src/atm_dyn_iconam/mo_solve_nonhydro.f90": [
        (465, 466),
        (497, 498)]}:
        raise NotImplementedError("Demo SDFG doesn't support this integration!")

    # create SDFG
    current_folder = os.path.dirname(os.path.abspath(__file__))

    demo_sdfg = dace.SDFG.from_file(current_folder + "/velocity_tendencies_autoopt_1.sdfgz")
    demo_sdfg.save(str(output_folder / "velocity_tendencies_unsimplified.sdfgz"))

    # create module definitions
    with open(output_folder / "velocity_tendencies_module_definitions.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            VELOCITY_TENDENCIES_MODULE_DEFINITIONS,
            Dumper=YAML_Dumper,
        ))

    # create associations
    with open(output_folder / "velocity_tendencies_associations.yaml", "w") as module_definitions_yaml:
        module_definitions_yaml.write(dump_yaml(
            VELOCITY_TENDENCIES_ASSOCIATIONS,
            Dumper=YAML_Dumper,
        ))



###############################################################################
# End demo SDFG
###############################################################################

def main():

    parser = argparse.ArgumentParser(description=(
        "Generate SDFGs and corresponding helper files (e.g., for integration)"
    ))

    parser.add_argument(
        "sdfg_name",
        help="The name of the SDFG to generate",
    )
    parser.add_argument(
        "integrations_yaml_path",
        help="Path to `integrations.yaml` file",
        type=Path,
    )
    parser.add_argument(
        "output_folder",
        help="Output folder for SDFG and helper files",
        type=Path,
    )
    args = parser.parse_args()

    with open(args.integrations_yaml_path) as yaml_file:
        integrations_yaml = load_yaml(yaml_file, Loader=YAML_Loader)

    required_associations = defaultdict(lambda: [])
    for fortran_source_file, fortran_source_integrations in integrations_yaml.items():
        fortran_source_file = str(Path(fortran_source_file))
        for sdfg_name, line_nr_integrations in fortran_source_integrations.items():
            if sdfg_name == args.sdfg_name:
                for line_nr_integration in line_nr_integrations:
                    start_line_nr, end_line_nr = line_nr_integration
                    required_associations[fortran_source_file].append((start_line_nr, end_line_nr))

    raise Exception(args.sdfg_name)
    if args.sdfg_name == "get_albedos":
        create_get_albedos_sdfg(args.output_folder, required_associations)
    elif args.sdfg_name == "edges_to_cells_bilinear_interpolation":
        create_edges_to_cells_bilinear_interpolation_sdfg(args.output_folder, required_associations)
    elif args.sdfg_name == "velocity_tendencies":
        create_velocity_tendencies_sdfg(args.output_folder, required_associations)
    else:
        parser.error(f"Unknown SDFG '{args.sdfg_name}'!")


if __name__ == "__main__":
    main()
