#include <dace/dace.h>

struct thermodynamics_type {

};
struct config_type {
    int __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8;
    int __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9;
    int __f2dace_SA_i_emiss_from_band_lw_d_0_s_5;
    int __f2dace_SA_sw_albedo_weights_d_0_s_6;
    int __f2dace_SA_sw_albedo_weights_d_1_s_7;
    int __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8;
    int __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9;
    int __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5;
    int __f2dace_SOA_sw_albedo_weights_d_0_s_6;
    int __f2dace_SOA_sw_albedo_weights_d_1_s_7;
    int* i_band_from_reordered_g_lw;
    int* i_band_from_reordered_g_sw;
    int* i_emiss_from_band_lw;
    double* sw_albedo_weights;
};
struct cloud_type {

};
struct flux_type {

};
struct aerosol_type {

};
struct gas_type {

};
struct single_level_type {
    int __f2dace_SA_lw_emissivity_d_0_s_20;
    int __f2dace_SA_lw_emissivity_d_1_s_21;
    int __f2dace_SA_sw_albedo_d_0_s_16;
    int __f2dace_SA_sw_albedo_d_1_s_17;
    int __f2dace_SA_sw_albedo_direct_d_0_s_18;
    int __f2dace_SA_sw_albedo_direct_d_1_s_19;
    int __f2dace_SOA_lw_emissivity_d_0_s_20;
    int __f2dace_SOA_lw_emissivity_d_1_s_21;
    int __f2dace_SOA_sw_albedo_d_0_s_16;
    int __f2dace_SOA_sw_albedo_d_1_s_17;
    int __f2dace_SOA_sw_albedo_direct_d_0_s_18;
    int __f2dace_SOA_sw_albedo_direct_d_1_s_19;
    double* lw_emissivity;
    double* sw_albedo;
    double* sw_albedo_direct;
};

typedef void * get_albedosHandle_t;
extern "C" get_albedosHandle_t __dace_init_get_albedos(int __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0, int __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0, int __f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0, int __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1, int __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1, int __f2dace_SA_sw_albedo_d_0_s_16_single_level_1, int __f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1, int __f2dace_SA_sw_albedo_weights_d_0_s_6_config_0, int __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0, int __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0, int __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0, int __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1, int __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1, int __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1, int __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1, int __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1, int __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1, int __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0, int __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0, int sym_iendcol, int sym_istartcol);
extern "C" int __dace_exit_get_albedos(get_albedosHandle_t handle);
extern "C" void __program_get_albedos(get_albedosHandle_t handle, aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, double * __restrict__ lw_albedo_var_599, single_level_type* single_level, double * __restrict__ sw_albedo_diffuse_var_601, double * __restrict__ sw_albedo_direct_var_600, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol);
