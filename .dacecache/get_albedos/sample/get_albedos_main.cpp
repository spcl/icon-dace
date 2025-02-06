#include <cstdlib>
#include "../include/get_albedos.h"

int main(int argc, char **argv) {
    get_albedosHandle_t handle;
    int err;
    int iendcol = 42;
    int istartcol = 42;
    int ncol = 42;
    int nlev = 42;
    int nulout = 42;
    int sym_iendcol = 42;
    int sym_istartcol = 42;
    double * __restrict__ lw_albedo_var_599 = (double*) calloc((((140 * sym_iendcol) - (140 * sym_istartcol)) + 140), sizeof(double));
    double * __restrict__ sw_albedo_diffuse_var_601 = (double*) calloc((((112 * sym_iendcol) - (112 * sym_istartcol)) + 112), sizeof(double));
    double * __restrict__ sw_albedo_direct_var_600 = (double*) calloc((((112 * sym_iendcol) - (112 * sym_istartcol)) + 112), sizeof(double));


    handle = __dace_init_get_albedos(__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0, __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0, __f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0, __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1, __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1, __f2dace_SA_sw_albedo_d_0_s_16_single_level_1, __f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1, __f2dace_SA_sw_albedo_weights_d_0_s_6_config_0, __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0, __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0, __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0, __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1, __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1, __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1, __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1, __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1, __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1, __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0, __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0, sym_iendcol, sym_istartcol);
    __program_get_albedos(handle, aerosol, cloud, config, flux, gas, lw_albedo_var_599, single_level, sw_albedo_diffuse_var_601, sw_albedo_direct_var_600, thermodynamics, iendcol, istartcol, ncol, nlev, nulout, sym_iendcol, sym_istartcol);
    err = __dace_exit_get_albedos(handle);

    free(lw_albedo_var_599);
    free(sw_albedo_diffuse_var_601);
    free(sw_albedo_direct_var_600);


    return err;
}
