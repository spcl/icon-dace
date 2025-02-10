#!/usr/bin/env python
import os
from typing import Dict, List, Tuple
from collections import defaultdict
from pathlib import Path
import argparse

from yaml import dump as dump_yaml, load as load_yaml
import yaml

try:
    from yaml import CDumper as YAML_Dumper, CLoader as YAML_Loader
except ImportError:
    from yaml import Dumper as YAML_Dumper, Loader as YAML_Loader

import dace


def create_radiation_sdfg(
    output_folder: Path,
    associations: Dict[str, List[Tuple[int, int]]],
):
    # create SDFG
    current_folder = os.path.dirname(os.path.abspath(__file__))
    demo_sdfg = dace.SDFG.from_file(current_folder + "/radiation.sdfgz")

    # sw_albedo_diffuse_var_601="0.0d0",
    # sw_albedo_direct_var_600="0.0d0",
    # lw_albedo_var_599="0.0d0",
    # Hotfix to compile
    for arr_name, arr in demo_sdfg.arrays.items():
        if (
            arr_name == "sw_albedo_diffuse_var_601"
            or arr_name == "sw_albedo_direct_var_600"
            or arr_name == "lw_albedo_var_599"
        ):
            assert arr.transient is False
            arr.transient = True

    demo_sdfg.save(str(output_folder / "radiation_unsimplified.sdfgz"))

    yaml_output = yaml.dump(RADIATION_MODULE_DEFINITIONS, default_flow_style=False)

    with open(output_folder / f"radiation_module_definitions.yam'", "w") as f:
        f.write(yaml_output)

    # create module definitions
    with open(
        output_folder / "radiation_module_definitions.yaml", "w"
    ) as module_definitions_yaml:
        module_definitions_yaml.write(
            dump_yaml(
                RADIATION_MODULE_DEFINITIONS,
                Dumper=YAML_Dumper,
            )
        )

    # create associations
    with open(
        output_folder / "radiation_associations.yaml", "w"
    ) as module_definitions_yaml:
        module_definitions_yaml.write(
            dump_yaml(
                RADIATION_ASSOCIATIONS,
                Dumper=YAML_Dumper,
            )
        )


VELOCITY_TENDENCIES_MODULE_DEFINITIONS = dict(
    # The only symbol needed for the actual SDFG
    t_patch="mo_model_domain",
    t_grid_cells="mo_model_domain",
    t_grid_edges="mo_model_domain",
    t_grid_vertices="mo_model_domain",
    t_subset_range="mo_model_domain",
    t_tangent_vectors="mo_model_domain",
    t_grid_geometry_info="mo_grid_geometry_info",
    t_geographical_coordinates="mo_math_types",
    t_cartesian_coordinates="mo_math_types",
    t_comm_pattern="mo_communication",
    t_comm_gather_pattern="mo_communication",
    t_scatterPattern="mo_communication",
    t_comm_pattern_collection="mo_communication",
    t_p_comm_pattern="mo_communication",
    t_uuid="mo_util_uuid_types",
    t_grid_domain_decomp_info="mo_decomposition_tools",
    t_glb2loc_index_lookup="mo_decomposition_tools",
    t_distrib_read_data="mo_read_netcdf_types",
    t_dist_dir="mo_dist_dir",
    t_int_state="mo_intp_data_strc",
    t_lsq="mo_intp_data_strc",
    t_cell_environ="mo_intp_data_strc",
    t_gauss_quad="mo_intp_data_strc",
    t_nh_state="mo_nonhydro_types",
    t_nh_diag="mo_nonhydro_types",
    t_nh_metrics="mo_nonhydro_types",
    t_nh_prog="mo_nonhydro_types",
    t_nh_ref="mo_nonhydro_types",
    t_prepare_adv="mo_prepadv_types",
    # helpers for fortan interface
    message="mo_exception",
    em_error="mo_exception",
    MAX_CHAR_LENGTH="mo_impl_constants",
    # globals
    dp="mo_kind",
    wp="mo_kind",
    sp="mo_kind",
    vp="mo_kind",
    i4="mo_kind",
    em_none="mo_exception",
    em_info="mo_exception",
    em_warn="mo_exception",
    em_param="mo_exception",
    em_debug="mo_exception",
    number_of_warnings="mo_exception",
    number_of_errors="mo_exception",
    warning="mo_exception",
    grid_length_rescale_factor="mo_grid_config",
    grid_rescale_factor="mo_grid_config",
    grid_angular_velocity="mo_grid_config",
    grid_sphere_radius="mo_grid_config",
    lrescale_ang_vel="mo_grid_config",
    namelist_grid_angular_velocity="mo_grid_config",
    no_of_dynamics_grids="mo_grid_config",
    n_dom="mo_grid_config",
    n_dom_start="mo_grid_config",
    sphere_geometry="mo_grid_geometry_info",
    triangular_cell="mo_grid_geometry_info",
    dbl_eps="mo_math_constants",
    rd="mo_physical_constants",
    cpd="mo_physical_constants",
    cvd="mo_physical_constants",
    cvd_o_rd="mo_physical_constants",
    earth_radius="mo_physical_constants",
    earth_angular_velocity="mo_physical_constants",
    alf="mo_physical_constants",
    als="mo_physical_constants",
    alsdcp="mo_physical_constants",
    alv="mo_physical_constants",
    alvdcp="mo_physical_constants",
    clw="mo_physical_constants",
    con0_h="mo_physical_constants",
    con_h="mo_physical_constants",
    con_m="mo_physical_constants",
    cpv="mo_physical_constants",
    cv_i="mo_physical_constants",
    cv_v="mo_physical_constants",
    cvv="mo_physical_constants",
    dv0="mo_physical_constants",
    eta0d="mo_physical_constants",
    grav="mo_physical_constants",
    o_m_rdv="mo_physical_constants",
    p0ref="mo_physical_constants",
    rcpd="mo_physical_constants",
    rcpl="mo_physical_constants",
    rcpv="mo_physical_constants",
    rcvd="mo_physical_constants",
    rd_o_cpd="mo_physical_constants",
    rdv="mo_physical_constants",
    rhoh2o="mo_physical_constants",
    rhoice="mo_physical_constants",
    rv="mo_physical_constants",
    t3="mo_physical_constants",
    tmelt="mo_physical_constants",
    vtmpc1="mo_physical_constants",
    vtmpc2="mo_physical_constants",
    filename_max="mo_io_units",
    nerr="mo_io_units",
    nlog="mo_io_units",
    mpi_comm_null="mo_mpi",
    p_mpi_comm_null="mo_mpi",
    comm_lev="mo_mpi",
    p_io="mo_mpi",
    p_pe="mo_mpi",
    proc_split="mo_mpi",
    mpi_request_null="mo_mpi",
    process_mpi_stdio_id="mo_mpi",
    i_am_accel_node="mo_mpi",
    max_hw="mo_impl_constants",
    min_rlcell_int="mo_impl_constants",
    min_rlcell="mo_impl_constants",
    min_rlvert_int="mo_impl_constants",
    min_rlvert="mo_impl_constants",
    min_rledge_int="mo_impl_constants",
    min_rledge="mo_impl_constants",
    max_dom="mo_impl_constants",
    rayleigh_classic="mo_impl_constants",
    rayleigh_klemp="mo_impl_constants",
    grf_bdywidth_c="mo_impl_constants_grf",
    grf_bdywidth_e="mo_impl_constants_grf",
    grf_bdywidth_v="mo_impl_constants_grf",
    nproma="mo_parallel_config",
    nblocks_c="mo_parallel_config",
    ignore_nproma_use_nblocks_c="mo_parallel_config",
    p_test_run="mo_parallel_config",
    l_test_openmp="mo_parallel_config",
    icon_comm_openmp="mo_parallel_config",
    num_io_procs="mo_parallel_config",
    num_restart_procs="mo_parallel_config",
    num_prefetch_proc="mo_parallel_config",
    proc0_shift="mo_parallel_config",
    use_dycore_barrier="mo_parallel_config",
    is_iau_active="mo_initicon_config",
    iau_wgt_dyn="mo_initicon_config",
    edge2cell_coeff_cc="mo_intp_data_strc",
)

sub_dict = dict(
    p_prog="p_nh%prog(nnew)",
    p_patch="p_patch",
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
    p_patch="p_patch",
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

RADIATION_MODULE_DEFINITIONS = {
    "ecrad": "mo_ecrad",
    "ecrad_ssi_default": "mo_ecrad",
    "ISolverSpartacus": "mo_ecrad",
    "t_ecrad_conf": "mo_ecrad",
    "t_ecrad_aerosol_type": "mo_ecrad",
    "t_ecrad_single_level_type": "mo_ecrad",
    "t_ecrad_thermodynamics_type": "mo_ecrad",
    "t_ecrad_gas_type": "mo_ecrad",
    "t_ecrad_flux_type": "mo_ecrad",
    "t_ecrad_cloud_type": "mo_ecrad",
    "config_type": "radiation_config",
    "aerosol_type": "radiation_aerosol",
    "single_level_type": "radiation_single_level",
    "thermodynamics_type": "radiation_thermodynamics",
    "gas_type": "radiation_gas",
    "flux_type": "radiation_flux",
    "cloud_type": "radiation_cloud",
    "t_opt_ptrs": "mo_ecrad",
    "ecrad_hyd_list": "mo_ecrad",
    "ecrad_iqr": "mo_ecrad",
    "ecrad_iqs": "mo_ecrad",
    "ecrad_iqg": "mo_ecrad",
    "nulout": "radiation_io",
    "ncol": "mo_parallel_config",
    "nproma_sub": "mo_parallel_config",
    "c_null_ptr": "iso_c_binding",
}

RADIATION_MODULE_DEFINITIONS = (
    VELOCITY_TENDENCIES_MODULE_DEFINITIONS | RADIATION_MODULE_DEFINITIONS
)

sub_dict3 = dict(
    aerosol="ecrad_aerosol",
    cloud="ecrad_cloud",
    config="ecrad_conf",
    flux="ecrad_flux",
    gas="ecrad_gas",
    single_level="ecrad_single_level",
    thermodynamics="ecrad_thermodynamics",
    iendcol="i_endidx_rad",
    istartcol="i_startidx_rad",
    ncol="nproma_sub",
    nlev="nlev",
    #sym_iendcol="i_endidx_rad",
    #sym_istartcol="i_startidx_rad",
    # sw_albedo_diffuse_var_601="0.0d0",
    # sw_albedo_direct_var_600="0.0d0",
    # lw_albedo_var_599="0.0d0",
)

VELOCITY_TENDENCIES_ASSOCIATIONS = {
    "mo_solve_nonhydro.f90": {465: {466: sub_dict}, 497: {498: sub_dict2}}
}

RADIATION_ASSOCIATIONS = {
    "mo_nwp_ecrad_interface.f90": {
        432: {442: sub_dict3},
    }
}


def create_velocity_tendencies_sdfg(
    output_folder: Path,
    associations: Dict[str, List[Tuple[int, int]]],
):

    if associations != {
        "src/atm_dyn_iconam/mo_solve_nonhydro.f90": [(465, 466), (497, 498)]
    }:
        raise NotImplementedError("Demo SDFG doesn't support this integration!")

    # create SDFG
    current_folder = os.path.dirname(os.path.abspath(__file__))

    demo_sdfg = dace.SDFG.from_file(
        current_folder + "/velocity_tendencies_autoopt_1.sdfgz"
    )
    demo_sdfg.save(str(output_folder / "velocity_tendencies_unsimplified.sdfgz"))

    # create module definitions
    with open(
        output_folder / "velocity_tendencies_module_definitions.yaml", "w"
    ) as module_definitions_yaml:
        module_definitions_yaml.write(
            dump_yaml(
                VELOCITY_TENDENCIES_MODULE_DEFINITIONS,
                Dumper=YAML_Dumper,
            )
        )

    # create associations
    with open(
        output_folder / "velocity_tendencies_associations.yaml", "w"
    ) as module_definitions_yaml:
        module_definitions_yaml.write(
            dump_yaml(
                VELOCITY_TENDENCIES_ASSOCIATIONS,
                Dumper=YAML_Dumper,
            )
        )


###############################################################################
# End demo SDFG
###############################################################################


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate SDFGs and corresponding helper files (e.g., for integration)"
        )
    )

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
                    required_associations[fortran_source_file].append(
                        (start_line_nr, end_line_nr)
                    )

    if args.sdfg_name == "radiation":
        create_radiation_sdfg(args.output_folder, required_associations)
    elif args.sdfg_name == "velocity_tendencies":
        create_velocity_tendencies_sdfg(args.output_folder, required_associations)
    else:
        parser.error(f"Unknown SDFG '{args.sdfg_name}'!")


if __name__ == "__main__":
    main()
