/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"
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
struct cloud_type {

};
struct aerosol_type {

};
struct flux_type {

};
struct gas_type {

};
struct thermodynamics_type {

};

struct radiation_state_t {

};

int __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0;
int __f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0;
int __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0;
int __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0;
int __f2dace_SA_sw_albedo_weights_d_0_s_6_config_0;
int __f2dace_SA_sw_albedo_weights_d_1_s_7_config_0;
int __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0;
int __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0;
int __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0;
int __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0;
int __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1;
int __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1;
int __f2dace_SA_sw_albedo_d_0_s_16_single_level_1;
int __f2dace_SA_sw_albedo_d_1_s_17_single_level_1;
int __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1;
int __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1;
int __f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1;
int __f2dace_SA_sw_albedo_direct_d_1_s_19_single_level_1;
int __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1;
int __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1;
int __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1;
int __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1;
void __program_radiation_internal(radiation_state_t*__state, aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, single_level_type* single_level, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    config_type** config_var_379_0;
    config_var_379_0 = &config;
    double* v_config_var_379_sw_albedo_weights;
    v_config_var_379_sw_albedo_weights = (double*)(&((*config_var_379_0)->sw_albedo_weights)[0]);
    int* v_config_var_379_i_band_from_reordered_g_sw;
    v_config_var_379_i_band_from_reordered_g_sw = (int*)(&((*config_var_379_0)->i_band_from_reordered_g_sw)[0]);
    int* v_config_var_379_i_emiss_from_band_lw;
    v_config_var_379_i_emiss_from_band_lw = (int*)(&((*config_var_379_0)->i_emiss_from_band_lw)[0]);
    int* v_config_var_379_i_band_from_reordered_g_lw;
    v_config_var_379_i_band_from_reordered_g_lw = (int*)(&((*config_var_379_0)->i_band_from_reordered_g_lw)[0]);
    double *sw_albedo_band_0 = nullptr;
    double *tmp_call_3_0;
    tmp_call_3_0 = new double DACE_ALIGN(64)[((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * (__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 - 1)) + __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1)];
    double *tmp_libnode_0_0;
    tmp_libnode_0_0 = new double DACE_ALIGN(64)[((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * (__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 - 1)) + __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1)];
    double *tmp_noncontig_0_0 = nullptr;
    double *tmp_noncontig_1_0 = nullptr;
    double *tmp_noncontig_2_0;
    tmp_noncontig_2_0 = new double DACE_ALIGN(64)[__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0];
    double *tmp_noncontig_3_0 = nullptr;
    single_level_type** this_var_378_0;
    this_var_378_0 = &single_level;
    int sym_istartcol_var_380_0;
    int sym_iendcol_var_381_0;
    int tmp_call_1_0;
    int tmp_index_38_0;
    int tmp_parfor_2_0_0;
    int tmp_parfor_1_1_0;
    int tmp_index_32_0;
    int tmp_index_29_0;
    int tmp_index_33_0;
    int tmp_parfor_1_0;
    int _for_it_2_0;
    int _if_cond_1_0;
    int _for_it_3_0;
    int tmp_index_6_0;
    int tmp_index_4_0;
    int _for_it_4_0;
    int tmp_parfor_0_0;
    int _for_it_0_0;
    int tmp_index_0_0;
    int _for_it_1_0;
    int _for_it_5_0;
    int tmp_index_17_0;
    int _for_it_6_0;
    int tmp_parfor_6_0;
    int tmp_index_49_0;
    int tmp_index_51_0;
    int tmp_parfor_5_0;
    int tmp_parfor_0_1_0;
    int tmp_index_12_0;
    int tmp_index_15_0;
    int tmp_index_16_0;
    int tmp_parfor_4_0;
    int tmp_parfor_3_0;
    int tmp_parfor_3_1_0;
    int tmp_index_39_0;
    int tmp_index_43_0;
    int tmp_parfor_2_1;
    int _for_it_7_0;
    int _if_cond_2_0;
    int _for_it_8_0;
    int tmp_index_23_0;
    int tmp_index_21_0;
    int _for_it_9_0;


    {

        {
            int nulout_out;

            ///////////////////
            // Tasklet code (T_l706_c706)
            nulout_out = 6;
            ///////////////////

            nulout = nulout_out;
        }

    }

    {

        {
            int nulout_out;

            ///////////////////
            // Tasklet code (T_l7809_c7809)
            nulout_out = 6;
            ///////////////////

            nulout = nulout_out;
        }

    }
    sym_istartcol_var_380_0 = istartcol;
    sym_iendcol_var_381_0 = iendcol;
    sw_albedo_band_0 = new double DACE_ALIGN(64)[(((14 * sym_iendcol_var_381_0) - (14 * sym_istartcol_var_380_0)) + 14)];

    for (_for_it_0_0 = 1; (_for_it_0_0 <= 14); _for_it_0_0 = (_for_it_0_0 + 1)) {
        for (_for_it_1_0 = istartcol; (_for_it_1_0 <= iendcol); _for_it_1_0 = (_for_it_1_0 + 1)) {

            tmp_index_0_0 = (_for_it_1_0 - istartcol);
            {

                {
                    double sw_albedo_band_out_0;

                    ///////////////////
                    // Tasklet code (T_l1003_c1003)
                    sw_albedo_band_out_0 = 0.0;
                    ///////////////////

                    sw_albedo_band_0[(tmp_index_0_0 + ((_for_it_0_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = sw_albedo_band_out_0;
                }

            }

        }

    }

    for (_for_it_2_0 = 1; (_for_it_2_0 <= 14); _for_it_2_0 = (_for_it_2_0 + 1)) {
        for (_for_it_3_0 = 1; (_for_it_3_0 <= 2); _for_it_3_0 = (_for_it_3_0 + 1)) {
            {


            }
            _if_cond_1_0 = (v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_2_0)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_3_0)] != 0.0);
            if ((_if_cond_1_0 == 1)) {
                for (_for_it_4_0 = istartcol; (_for_it_4_0 <= iendcol); _for_it_4_0 = (_for_it_4_0 + 1)) {

                    tmp_index_6_0 = (_for_it_4_0 - istartcol);
                    tmp_index_4_0 = (_for_it_4_0 - istartcol);
                    {
                        double* v_this_var_378_sw_albedo;
                        v_this_var_378_sw_albedo = (double*)(&((*this_var_378_0)->sw_albedo)[0]);

                        {
                            double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_2_0)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_3_0)];
                            double sw_albedo_band_0_in_0 = sw_albedo_band_0[(tmp_index_6_0 + ((_for_it_2_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                            double this_var_378_0_in_sw_albedo_0 = v_this_var_378_sw_albedo[(((__f2dace_SA_sw_albedo_d_0_s_16_single_level_1 * ((- __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1) + _for_it_3_0)) - __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1) + _for_it_4_0)];
                            double sw_albedo_band_out_0;

                            ///////////////////
                            // Tasklet code (T_l1010_c1010)
                            sw_albedo_band_out_0 = (sw_albedo_band_0_in_0 + (config_var_379_0_in_sw_albedo_weights_0 * this_var_378_0_in_sw_albedo_0));
                            ///////////////////

                            sw_albedo_band_0[(tmp_index_4_0 + ((_for_it_2_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = sw_albedo_band_out_0;
                        }

                    }

                }
            }

        }

    }
    tmp_noncontig_0_0 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + ((__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1))) + 1)];

    for (tmp_parfor_0_1_0 = 1; (tmp_parfor_0_1_0 <= __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0); tmp_parfor_0_1_0 = (tmp_parfor_0_1_0 + 1)) {
        for (tmp_parfor_0_0 = istartcol; (tmp_parfor_0_0 <= iendcol); tmp_parfor_0_0 = (tmp_parfor_0_0 + 1)) {

            tmp_index_12_0 = (tmp_parfor_0_0 - istartcol);
            tmp_index_15_0 = ((tmp_parfor_0_0 + (istartcol - istartcol)) - istartcol);
            {


            }
            tmp_index_16_0 = (v_config_var_379_i_band_from_reordered_g_sw[((- __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0) + tmp_parfor_0_1_0)] - 1);
            {

                {
                    double sw_albedo_band_0_in_0 = sw_albedo_band_0[(tmp_index_15_0 + (tmp_index_16_0 * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                    double tmp_noncontig_0_out_0;

                    ///////////////////
                    // Tasklet code (T_l1015_c1015)
                    tmp_noncontig_0_out_0 = sw_albedo_band_0_in_0;
                    ///////////////////

                    tmp_noncontig_0_0[(tmp_index_12_0 + ((tmp_parfor_0_1_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = tmp_noncontig_0_out_0;
                }

            }

        }

    }

    {
        double *sw_albedo_diffuse_var_601;
        sw_albedo_diffuse_var_601 = new double DACE_ALIGN(64)[(((112 * sym_iendcol) - (112 * sym_istartcol)) + 112)];

        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1); __i0 += 1) {
                for (auto __i1 = 0; __i1 < __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0; __i1 += 1) {
                    {
                        double __inp = tmp_noncontig_0_0[(__i0 + (__i1 * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                        double __out;

                        ///////////////////
                        // Tasklet code (transpose)
                        __out = __inp;
                        ///////////////////

                        sw_albedo_diffuse_var_601[((112 * __i0) + __i1)] = __out;
                    }
                }
            }
        }
        delete[] sw_albedo_diffuse_var_601;
        delete[] tmp_noncontig_0_0;

    }

    for (_for_it_5_0 = 1; (_for_it_5_0 <= 14); _for_it_5_0 = (_for_it_5_0 + 1)) {
        for (_for_it_6_0 = istartcol; (_for_it_6_0 <= iendcol); _for_it_6_0 = (_for_it_6_0 + 1)) {

            tmp_index_17_0 = (_for_it_6_0 - istartcol);
            {

                {
                    double sw_albedo_band_out_0;

                    ///////////////////
                    // Tasklet code (T_l1018_c1018)
                    sw_albedo_band_out_0 = 0.0;
                    ///////////////////

                    sw_albedo_band_0[(tmp_index_17_0 + ((_for_it_5_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = sw_albedo_band_out_0;
                }

            }

        }

    }

    for (_for_it_7_0 = 1; (_for_it_7_0 <= 14); _for_it_7_0 = (_for_it_7_0 + 1)) {
        for (_for_it_8_0 = 1; (_for_it_8_0 <= 2); _for_it_8_0 = (_for_it_8_0 + 1)) {
            {


            }
            _if_cond_2_0 = (v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_7_0)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_8_0)] != 0.0);
            if ((_if_cond_2_0 == 1)) {
                for (_for_it_9_0 = istartcol; (_for_it_9_0 <= iendcol); _for_it_9_0 = (_for_it_9_0 + 1)) {

                    tmp_index_23_0 = (_for_it_9_0 - istartcol);
                    tmp_index_21_0 = (_for_it_9_0 - istartcol);
                    {
                        double* v_this_var_378_sw_albedo_direct;
                        v_this_var_378_sw_albedo_direct = (double*)(&((*this_var_378_0)->sw_albedo_direct)[0]);

                        {
                            double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_7_0)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_8_0)];
                            double sw_albedo_band_0_in_0 = sw_albedo_band_0[(tmp_index_23_0 + ((_for_it_7_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                            double this_var_378_0_in_sw_albedo_direct_0 = v_this_var_378_sw_albedo_direct[(((__f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1 * ((- __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1) + _for_it_8_0)) - __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1) + _for_it_9_0)];
                            double sw_albedo_band_out_0;

                            ///////////////////
                            // Tasklet code (T_l1025_c1025)
                            sw_albedo_band_out_0 = (sw_albedo_band_0_in_0 + (config_var_379_0_in_sw_albedo_weights_0 * this_var_378_0_in_sw_albedo_direct_0));
                            ///////////////////

                            sw_albedo_band_0[(tmp_index_21_0 + ((_for_it_7_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = sw_albedo_band_out_0;
                        }

                    }

                }
            }

        }

    }
    tmp_noncontig_1_0 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + ((__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1))) + 1)];

    for (tmp_parfor_1_1_0 = 1; (tmp_parfor_1_1_0 <= __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0); tmp_parfor_1_1_0 = (tmp_parfor_1_1_0 + 1)) {
        for (tmp_parfor_1_0 = istartcol; (tmp_parfor_1_0 <= iendcol); tmp_parfor_1_0 = (tmp_parfor_1_0 + 1)) {

            tmp_index_32_0 = ((tmp_parfor_1_0 + (istartcol - istartcol)) - istartcol);
            tmp_index_29_0 = (tmp_parfor_1_0 - istartcol);
            {


            }
            tmp_index_33_0 = (v_config_var_379_i_band_from_reordered_g_sw[((- __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0) + tmp_parfor_1_1_0)] - 1);
            {

                {
                    double sw_albedo_band_0_in_0 = sw_albedo_band_0[(tmp_index_32_0 + (tmp_index_33_0 * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                    double tmp_noncontig_1_out_0;

                    ///////////////////
                    // Tasklet code (T_l1030_c1030)
                    tmp_noncontig_1_out_0 = sw_albedo_band_0_in_0;
                    ///////////////////

                    tmp_noncontig_1_0[(tmp_index_29_0 + ((tmp_parfor_1_1_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = tmp_noncontig_1_out_0;
                }

            }

        }

    }
    tmp_call_1_0 = -2147483648LL;

    delete[] sw_albedo_band_0;
    {
        double *sw_albedo_direct_var_600;
        sw_albedo_direct_var_600 = new double DACE_ALIGN(64)[(((112 * sym_iendcol) - (112 * sym_istartcol)) + 112)];

        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1); __i0 += 1) {
                for (auto __i1 = 0; __i1 < __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0; __i1 += 1) {
                    {
                        double __inp = tmp_noncontig_1_0[(__i0 + (__i1 * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                        double __out;

                        ///////////////////
                        // Tasklet code (transpose)
                        __out = __inp;
                        ///////////////////

                        sw_albedo_direct_var_600[((112 * __i0) + __i1)] = __out;
                    }
                }
            }
        }
        delete[] sw_albedo_direct_var_600;
        delete[] tmp_noncontig_1_0;

    }

    for (tmp_parfor_0_0 = 1; (tmp_parfor_0_0 <= ((__f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0 + 1) - 1)); tmp_parfor_0_0 = (tmp_parfor_0_0 + 1)) {
        {


        }
        if ((v_config_var_379_i_emiss_from_band_lw[((- __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0) + tmp_parfor_0_0)] > tmp_call_1_0)) {
            {


            }
            tmp_call_1_0 = v_config_var_379_i_emiss_from_band_lw[((- __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0) + tmp_parfor_0_0)];

        }

    }

    for (tmp_parfor_2_0_0 = 1; (tmp_parfor_2_0_0 <= __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0); tmp_parfor_2_0_0 = (tmp_parfor_2_0_0 + 1)) {
        {


        }
        tmp_index_38_0 = (v_config_var_379_i_band_from_reordered_g_lw[((- __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0) + tmp_parfor_2_0_0)] - __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0);
        {

            {
                int config_var_379_0_in_i_emiss_from_band_lw_0 = v_config_var_379_i_emiss_from_band_lw[tmp_index_38_0];
                double tmp_noncontig_2_out_0;

                ///////////////////
                // Tasklet code (T_l1035_c1035)
                tmp_noncontig_2_out_0 = config_var_379_0_in_i_emiss_from_band_lw_0;
                ///////////////////

                tmp_noncontig_2_0[(tmp_parfor_2_0_0 - 1)] = tmp_noncontig_2_out_0;
            }

        }

    }
    tmp_noncontig_3_0 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + ((__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1))) + 1)];

    for (tmp_parfor_3_1_0 = 1; (tmp_parfor_3_1_0 <= __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0); tmp_parfor_3_1_0 = (tmp_parfor_3_1_0 + 1)) {
        for (tmp_parfor_2_1 = istartcol; (tmp_parfor_2_1 <= iendcol); tmp_parfor_2_1 = (tmp_parfor_2_1 + 1)) {

            tmp_index_39_0 = (tmp_parfor_2_1 - istartcol);
            tmp_index_43_0 = (tmp_noncontig_2_0[(tmp_parfor_3_1_0 - 1)] - __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1);
            {
                double* v_this_var_378_lw_emissivity;
                v_this_var_378_lw_emissivity = (double*)(&((*this_var_378_0)->lw_emissivity)[0]);

                {
                    double this_var_378_0_in_lw_emissivity_0 = v_this_var_378_lw_emissivity[(((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * tmp_index_43_0) - __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_2_1)];
                    double tmp_noncontig_3_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    tmp_noncontig_3_out_0 = this_var_378_0_in_lw_emissivity_0;
                    ///////////////////

                    tmp_noncontig_3_0[(tmp_index_39_0 + ((tmp_parfor_3_1_0 - 1) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))] = tmp_noncontig_3_out_0;
                }

            }

        }

    }

    for (tmp_parfor_4_0 = __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1; (tmp_parfor_4_0 <= ((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 + __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) - 1)); tmp_parfor_4_0 = (tmp_parfor_4_0 + 1)) {
        for (tmp_parfor_3_0 = __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1; (tmp_parfor_3_0 <= ((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 + __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) - 1)); tmp_parfor_3_0 = (tmp_parfor_3_0 + 1)) {
            {

                {
                    double tmp_noncontig_3_0_in_0 = tmp_noncontig_3_0[(((- __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_3_0) + (((- __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) + tmp_parfor_4_0) * ((sym_iendcol_var_381_0 - sym_istartcol_var_380_0) + 1)))];
                    double tmp_libnode_0_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    tmp_libnode_0_out_0 = tmp_noncontig_3_0_in_0;
                    ///////////////////

                    tmp_libnode_0_0[(((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * ((- __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) + tmp_parfor_4_0)) - __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_3_0)] = tmp_libnode_0_out_0;
                }

            }

        }

    }

    delete[] tmp_noncontig_3_0;
    {

        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1; __i0 += 1) {
                for (auto __i1 = 0; __i1 < __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1; __i1 += 1) {
                    {
                        double __inp = tmp_libnode_0_0[((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * __i1) + __i0)];
                        double __out;

                        ///////////////////
                        // Tasklet code (transpose)
                        __out = __inp;
                        ///////////////////

                        tmp_call_3_0[((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * __i0) + __i1)] = __out;
                    }
                }
            }
        }

    }

    for (tmp_parfor_6_0 = istartcol; (tmp_parfor_6_0 <= ((((sym_iendcol - sym_istartcol) + 1) + istartcol) - 1)); tmp_parfor_6_0 = (tmp_parfor_6_0 + 1)) {
        for (tmp_parfor_5_0 = 1; (tmp_parfor_5_0 <= 140); tmp_parfor_5_0 = (tmp_parfor_5_0 + 1)) {

            tmp_index_49_0 = (tmp_parfor_6_0 - istartcol);
            tmp_index_51_0 = ((tmp_parfor_6_0 + (1 - istartcol)) - 1);
            {
                double *lw_albedo_var_599;
                lw_albedo_var_599 = new double DACE_ALIGN(64)[(((140 * sym_iendcol) - (140 * sym_istartcol)) + 140)];

                {
                    double tmp_call_3_0_in_0 = tmp_call_3_0[(((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * tmp_index_51_0) + tmp_parfor_5_0) - 1)];
                    double lw_albedo_var_382_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    lw_albedo_var_382_out_0 = (1.0 - tmp_call_3_0_in_0);
                    ///////////////////

                    lw_albedo_var_599[(((140 * tmp_index_49_0) + tmp_parfor_5_0) - 1)] = lw_albedo_var_382_out_0;
                }
                delete[] lw_albedo_var_599;

            }

        }

    }

    delete[] tmp_call_3_0;
    delete[] tmp_libnode_0_0;
    delete[] tmp_noncontig_2_0;
}

DACE_EXPORTED void __program_radiation(radiation_state_t *__state, aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, single_level_type* single_level, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    __program_radiation_internal(__state, aerosol, cloud, config, flux, gas, single_level, thermodynamics, iendcol, istartcol, ncol, nlev, nulout, sym_iendcol, sym_istartcol);
}

DACE_EXPORTED radiation_state_t *__dace_init_radiation(aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, single_level_type* single_level, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    int __result = 0;
    radiation_state_t *__state = new radiation_state_t;


    __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0 = config->__f2dace_SOA_i_emiss_from_band_lw_d_0_s_5;
    __f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0 = config->__f2dace_SA_i_emiss_from_band_lw_d_0_s_5;
    __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0 = config->__f2dace_SOA_sw_albedo_weights_d_0_s_6;
    __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0 = config->__f2dace_SOA_sw_albedo_weights_d_1_s_7;
    __f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 = config->__f2dace_SA_sw_albedo_weights_d_0_s_6;
    __f2dace_SA_sw_albedo_weights_d_1_s_7_config_0 = config->__f2dace_SA_sw_albedo_weights_d_1_s_7;
    __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0 = config->__f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8;
    __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0 = config->__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8;
    __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0 = config->__f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9;
    __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0 = config->__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9;
    __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1 = single_level->__f2dace_SOA_sw_albedo_d_0_s_16;
    __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1 = single_level->__f2dace_SOA_sw_albedo_d_1_s_17;
    __f2dace_SA_sw_albedo_d_0_s_16_single_level_1 = single_level->__f2dace_SA_sw_albedo_d_0_s_16;
    __f2dace_SA_sw_albedo_d_1_s_17_single_level_1 = single_level->__f2dace_SA_sw_albedo_d_1_s_17;
    __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1 = single_level->__f2dace_SOA_sw_albedo_direct_d_0_s_18;
    __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1 = single_level->__f2dace_SOA_sw_albedo_direct_d_1_s_19;
    __f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1 = single_level->__f2dace_SA_sw_albedo_direct_d_0_s_18;
    __f2dace_SA_sw_albedo_direct_d_1_s_19_single_level_1 = single_level->__f2dace_SA_sw_albedo_direct_d_1_s_19;
    __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1 = single_level->__f2dace_SOA_lw_emissivity_d_0_s_20;
    __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1 = single_level->__f2dace_SOA_lw_emissivity_d_1_s_21;
    __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 = single_level->__f2dace_SA_lw_emissivity_d_0_s_20;
    __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 = single_level->__f2dace_SA_lw_emissivity_d_1_s_21;

    if (__result) {
        delete __state;
        return nullptr;
    }
    return __state;
}

DACE_EXPORTED int __dace_exit_radiation(radiation_state_t *__state)
{
    int __err = 0;
    delete __state;
    return __err;
}
