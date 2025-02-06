#include <cstdlib>
#include "../include/radiation.h"

int main(int argc, char **argv) {
    radiationHandle_t handle;
    int err;
    int iendcol = 42;
    int istartcol = 42;
    int ncol = 42;
    int nlev = 42;
    int nulout = 42;
    int sym_iendcol = 42;
    int sym_istartcol = 42;


    handle = __dace_init_radiation(__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0, __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0, __f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0, __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1, __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1, __f2dace_SA_sw_albedo_d_0_s_16_single_level_1, __f2dace_SA_sw_albedo_d_1_s_17_single_level_1, __f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1, __f2dace_SA_sw_albedo_direct_d_1_s_19_single_level_1, __f2dace_SA_sw_albedo_weights_d_0_s_6_config_0, __f2dace_SA_sw_albedo_weights_d_1_s_7_config_0, __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0, __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0, __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0, __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1, __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1, __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1, __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1, __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1, __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1, __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0, __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0, sym_iendcol, sym_istartcol);
    __program_radiation(handle, aerosol, cloud, config, flux, gas, single_level, thermodynamics, iendcol, istartcol, ncol, nlev, nulout, sym_iendcol, sym_istartcol);
    err = __dace_exit_radiation(handle);



    return err;
}
