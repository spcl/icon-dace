/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"
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

struct get_albedos_state_t {

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
inline void global_init_fn0_0_30_0(get_albedos_state_t *__state, int&  nulout) {


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

}

inline void transpose_sdfg_2_3_2(get_albedos_state_t *__state, double* __restrict__ _inp, double* __restrict__ _out, int sym_iendcol_var_381, int sym_istartcol_var_380) {


    {

        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1); __i0 += 1) {
                for (auto __i1 = 0; __i1 < __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0; __i1 += 1) {
                    {
                        double __inp = _inp[(__i0 + (__i1 * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                        double __out;

                        ///////////////////
                        // Tasklet code (transpose)
                        __out = __inp;
                        ///////////////////

                        _out[((112 * __i0) + __i1)] = __out;
                    }
                }
            }
        }

    }

}

inline void transpose_sdfg_2_12_2(get_albedos_state_t *__state, double* __restrict__ _inp, double* __restrict__ _out) {


    {

        {
            #pragma omp parallel for
            for (auto __i0 = 0; __i0 < __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1; __i0 += 1) {
                for (auto __i1 = 0; __i1 < __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1; __i1 += 1) {
                    {
                        double __inp = _inp[((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * __i1) + __i0)];
                        double __out;

                        ///////////////////
                        // Tasklet code (transpose)
                        __out = __inp;
                        ///////////////////

                        _out[((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * __i0) + __i1)] = __out;
                    }
                }
            }
        }

    }

}

inline void get_albedos1_0_31_0(get_albedos_state_t *__state, config_type* __restrict__ config_var_379, const int* __restrict__ iendcol_var_381, const int* __restrict__ istartcol_var_380, single_level_type* __restrict__ this_var_378, double* __restrict__ lw_albedo_var_382, double* __restrict__ sw_albedo_diffuse_var_384, double* __restrict__ sw_albedo_direct_var_383, int sym_iendcol, int sym_iendcol_var_381, int sym_istartcol, int sym_istartcol_var_380) {
    double *sw_albedo_band;
    sw_albedo_band = new double DACE_ALIGN(64)[(((14 * sym_iendcol_var_381) - (14 * sym_istartcol_var_380)) + 14)];
    int tmp_call_1;
    double *tmp_call_3;
    tmp_call_3 = new double DACE_ALIGN(64)[((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * (__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 - 1)) + __f2dace_SA_lw_emissivity_d_1_s_21_single_level_1)];
    double *tmp_libnode_0;
    tmp_libnode_0 = new double DACE_ALIGN(64)[((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * (__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 - 1)) + __f2dace_SA_lw_emissivity_d_0_s_20_single_level_1)];
    double *tmp_noncontig_0;
    tmp_noncontig_0 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381 - sym_istartcol_var_380) + ((__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1))) + 1)];
    double *tmp_noncontig_1;
    tmp_noncontig_1 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381 - sym_istartcol_var_380) + ((__f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1))) + 1)];
    double *tmp_noncontig_2;
    tmp_noncontig_2 = new double DACE_ALIGN(64)[__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0];
    double *tmp_noncontig_3;
    tmp_noncontig_3 = new double DACE_ALIGN(64)[(((sym_iendcol_var_381 - sym_istartcol_var_380) + ((__f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1))) + 1)];
    int _if_cond_1;
    int _if_cond_2;
    double* v_config_var_379_sw_albedo_weights;
    v_config_var_379_sw_albedo_weights = (double*)(&(config_var_379->sw_albedo_weights)[0]);
    int* v_config_var_379_i_band_from_reordered_g_sw;
    v_config_var_379_i_band_from_reordered_g_sw = (int*)(&(config_var_379->i_band_from_reordered_g_sw)[0]);
    int* v_config_var_379_i_emiss_from_band_lw;
    v_config_var_379_i_emiss_from_band_lw = (int*)(&(config_var_379->i_emiss_from_band_lw)[0]);
    int _for_it_0;
    int tmp_index_0;
    int _for_it_1;
    int _for_it_2;
    int _for_it_3;
    int tmp_index_4;
    int tmp_index_6;
    int _for_it_4;
    int tmp_parfor_0_1;
    int tmp_index_12;
    int tmp_index_15;
    int tmp_parfor_0;
    int _for_it_5;
    int tmp_index_17;
    int _for_it_6;
    int _for_it_7;
    int _for_it_8;
    int tmp_index_21;
    int tmp_index_23;
    int _for_it_9;
    int tmp_parfor_1_1;
    int tmp_index_29;
    int tmp_index_32;
    int tmp_parfor_1;
    int tmp_parfor_2_0;
    int tmp_parfor_3_1;
    int tmp_index_39;
    int tmp_index_43;
    int tmp_parfor_2;
    int tmp_parfor_4;
    int tmp_parfor_3;
    int tmp_parfor_6;
    int tmp_index_49;
    int tmp_index_51;
    int tmp_parfor_5;


    for (_for_it_0 = 1; (_for_it_0 <= 14); _for_it_0 = (_for_it_0 + 1)) {
        for (_for_it_1 = istartcol_var_380; (_for_it_1 <= iendcol_var_381); _for_it_1 = (_for_it_1 + 1)) {

            tmp_index_0 = (_for_it_1 - istartcol_var_380[0]);
            {

                {
                    double sw_albedo_band_out_0;

                    ///////////////////
                    // Tasklet code (T_l1003_c1003)
                    sw_albedo_band_out_0 = 0.0;
                    ///////////////////

                    sw_albedo_band[(tmp_index_0 + ((_for_it_0 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = sw_albedo_band_out_0;
                }

            }

        }

    }

    for (_for_it_2 = 1; (_for_it_2 <= 14); _for_it_2 = (_for_it_2 + 1)) {
        for (_for_it_3 = 1; (_for_it_3 <= 2); _for_it_3 = (_for_it_3 + 1)) {
            {

                {
                    double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_2)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_3)];
                    int _if_cond_1_out;

                    ///////////////////
                    // Tasklet code (T_l0_c0)
                    _if_cond_1_out = (config_var_379_0_in_sw_albedo_weights_0 != 0.0);
                    ///////////////////

                    _if_cond_1 = _if_cond_1_out;
                }

            }
            if ((_if_cond_1 == 1)) {
                for (_for_it_4 = istartcol_var_380; (_for_it_4 <= iendcol_var_381); _for_it_4 = (_for_it_4 + 1)) {

                    tmp_index_4 = (_for_it_4 - istartcol_var_380[0]);
                    tmp_index_6 = (_for_it_4 - istartcol_var_380[0]);
                    {
                        double* v_this_var_378_sw_albedo;
                        v_this_var_378_sw_albedo = (double*)(&(this_var_378->sw_albedo)[0]);

                        {
                            double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_2)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_3)];
                            double sw_albedo_band_0_in_0 = sw_albedo_band[(tmp_index_6 + ((_for_it_2 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                            double this_var_378_0_in_sw_albedo_0 = v_this_var_378_sw_albedo[(((__f2dace_SA_sw_albedo_d_0_s_16_single_level_1 * ((- __f2dace_SOA_sw_albedo_d_1_s_17_single_level_1) + _for_it_3)) - __f2dace_SOA_sw_albedo_d_0_s_16_single_level_1) + _for_it_4)];
                            double sw_albedo_band_out_0;

                            ///////////////////
                            // Tasklet code (T_l1010_c1010)
                            sw_albedo_band_out_0 = (sw_albedo_band_0_in_0 + (config_var_379_0_in_sw_albedo_weights_0 * this_var_378_0_in_sw_albedo_0));
                            ///////////////////

                            sw_albedo_band[(tmp_index_4 + ((_for_it_2 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = sw_albedo_band_out_0;
                        }

                    }

                }
            }

        }

    }

    for (tmp_parfor_0_1 = 1; (tmp_parfor_0_1 <= __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0); tmp_parfor_0_1 = (tmp_parfor_0_1 + 1)) {
        for (tmp_parfor_0 = istartcol_var_380; (tmp_parfor_0 <= iendcol_var_381); tmp_parfor_0 = (tmp_parfor_0 + 1)) {

            tmp_index_12 = (tmp_parfor_0 - istartcol_var_380[0]);
            tmp_index_15 = ((tmp_parfor_0 + (istartcol_var_380[0] - istartcol_var_380[0])) - istartcol_var_380[0]);
            {
                int tmp_index_16;

                {
                    int config_var_379_0_in_i_band_from_reordered_g_sw_0 = v_config_var_379_i_band_from_reordered_g_sw[((- __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0) + tmp_parfor_0_1)];
                    int tmp_index_16_out;

                    ///////////////////
                    // Tasklet code (T_l1015_c1015)
                    tmp_index_16_out = (config_var_379_0_in_i_band_from_reordered_g_sw_0 - 1);
                    ///////////////////

                    tmp_index_16 = tmp_index_16_out;
                }
                {
                    double* sw_albedo_band_0_in = &sw_albedo_band[0];
                    int tmp_index_16_0_in = tmp_index_16;
                    double tmp_noncontig_0_out_0;

                    ///////////////////
                    // Tasklet code (T_l1015_c1015)
                    tmp_noncontig_0_out_0 = sw_albedo_band_0_in[(tmp_index_15 + (tmp_index_16_0_in * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                    ///////////////////

                    tmp_noncontig_0[(tmp_index_12 + ((tmp_parfor_0_1 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = tmp_noncontig_0_out_0;
                }

            }

        }

    }

    {

        transpose_sdfg_2_3_2(__state, &tmp_noncontig_0[0], &sw_albedo_diffuse_var_384[0], sym_iendcol_var_381, sym_istartcol_var_380);

    }

    for (_for_it_5 = 1; (_for_it_5 <= 14); _for_it_5 = (_for_it_5 + 1)) {
        for (_for_it_6 = istartcol_var_380; (_for_it_6 <= iendcol_var_381); _for_it_6 = (_for_it_6 + 1)) {

            tmp_index_17 = (_for_it_6 - istartcol_var_380[0]);
            {

                {
                    double sw_albedo_band_out_0;

                    ///////////////////
                    // Tasklet code (T_l1018_c1018)
                    sw_albedo_band_out_0 = 0.0;
                    ///////////////////

                    sw_albedo_band[(tmp_index_17 + ((_for_it_5 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = sw_albedo_band_out_0;
                }

            }

        }

    }

    for (_for_it_7 = 1; (_for_it_7 <= 14); _for_it_7 = (_for_it_7 + 1)) {
        for (_for_it_8 = 1; (_for_it_8 <= 2); _for_it_8 = (_for_it_8 + 1)) {
            {

                {
                    double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_7)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_8)];
                    int _if_cond_2_out;

                    ///////////////////
                    // Tasklet code (T_l0_c0)
                    _if_cond_2_out = (config_var_379_0_in_sw_albedo_weights_0 != 0.0);
                    ///////////////////

                    _if_cond_2 = _if_cond_2_out;
                }

            }
            if ((_if_cond_2 == 1)) {
                for (_for_it_9 = istartcol_var_380; (_for_it_9 <= iendcol_var_381); _for_it_9 = (_for_it_9 + 1)) {

                    tmp_index_21 = (_for_it_9 - istartcol_var_380[0]);
                    tmp_index_23 = (_for_it_9 - istartcol_var_380[0]);
                    {
                        double* v_this_var_378_sw_albedo_direct;
                        v_this_var_378_sw_albedo_direct = (double*)(&(this_var_378->sw_albedo_direct)[0]);

                        {
                            double config_var_379_0_in_sw_albedo_weights_0 = v_config_var_379_sw_albedo_weights[(((__f2dace_SA_sw_albedo_weights_d_0_s_6_config_0 * ((- __f2dace_SOA_sw_albedo_weights_d_1_s_7_config_0) + _for_it_7)) - __f2dace_SOA_sw_albedo_weights_d_0_s_6_config_0) + _for_it_8)];
                            double sw_albedo_band_0_in_0 = sw_albedo_band[(tmp_index_23 + ((_for_it_7 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                            double this_var_378_0_in_sw_albedo_direct_0 = v_this_var_378_sw_albedo_direct[(((__f2dace_SA_sw_albedo_direct_d_0_s_18_single_level_1 * ((- __f2dace_SOA_sw_albedo_direct_d_1_s_19_single_level_1) + _for_it_8)) - __f2dace_SOA_sw_albedo_direct_d_0_s_18_single_level_1) + _for_it_9)];
                            double sw_albedo_band_out_0;

                            ///////////////////
                            // Tasklet code (T_l1025_c1025)
                            sw_albedo_band_out_0 = (sw_albedo_band_0_in_0 + (config_var_379_0_in_sw_albedo_weights_0 * this_var_378_0_in_sw_albedo_direct_0));
                            ///////////////////

                            sw_albedo_band[(tmp_index_21 + ((_for_it_7 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = sw_albedo_band_out_0;
                        }

                    }

                }
            }

        }

    }

    for (tmp_parfor_1_1 = 1; (tmp_parfor_1_1 <= __f2dace_SA_i_band_from_reordered_g_sw_d_0_s_9_config_0); tmp_parfor_1_1 = (tmp_parfor_1_1 + 1)) {
        for (tmp_parfor_1 = istartcol_var_380; (tmp_parfor_1 <= iendcol_var_381); tmp_parfor_1 = (tmp_parfor_1 + 1)) {

            tmp_index_29 = (tmp_parfor_1 - istartcol_var_380[0]);
            tmp_index_32 = ((tmp_parfor_1 + (istartcol_var_380[0] - istartcol_var_380[0])) - istartcol_var_380[0]);
            {
                int tmp_index_33;

                {
                    int config_var_379_0_in_i_band_from_reordered_g_sw_0 = v_config_var_379_i_band_from_reordered_g_sw[((- __f2dace_SOA_i_band_from_reordered_g_sw_d_0_s_9_config_0) + tmp_parfor_1_1)];
                    int tmp_index_33_out;

                    ///////////////////
                    // Tasklet code (T_l1030_c1030)
                    tmp_index_33_out = (config_var_379_0_in_i_band_from_reordered_g_sw_0 - 1);
                    ///////////////////

                    tmp_index_33 = tmp_index_33_out;
                }
                {
                    double* sw_albedo_band_0_in = &sw_albedo_band[0];
                    int tmp_index_33_0_in = tmp_index_33;
                    double tmp_noncontig_1_out_0;

                    ///////////////////
                    // Tasklet code (T_l1030_c1030)
                    tmp_noncontig_1_out_0 = sw_albedo_band_0_in[(tmp_index_32 + (tmp_index_33_0_in * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                    ///////////////////

                    tmp_noncontig_1[(tmp_index_29 + ((tmp_parfor_1_1 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = tmp_noncontig_1_out_0;
                }

            }

        }

    }

    {

        transpose_sdfg_2_3_2(__state, &tmp_noncontig_1[0], &sw_albedo_direct_var_383[0], sym_iendcol_var_381, sym_istartcol_var_380);
        {
            int tmp_call_1_out;

            ///////////////////
            // Tasklet code (T_l0_c0)
            tmp_call_1_out = -2147483648LL;
            ///////////////////

            tmp_call_1 = tmp_call_1_out;
        }

    }

    for (tmp_parfor_0 = 1; (tmp_parfor_0 <= ((__f2dace_SA_i_emiss_from_band_lw_d_0_s_5_config_0 + 1) - 1)); tmp_parfor_0 = (tmp_parfor_0 + 1)) {
        {


        }
        if ((v_config_var_379_i_emiss_from_band_lw[((- __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0) + tmp_parfor_0)] > tmp_call_1)) {
            {

                {
                    int config_var_379_0_in_i_emiss_from_band_lw_0 = v_config_var_379_i_emiss_from_band_lw[((- __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0) + tmp_parfor_0)];
                    int tmp_call_1_out;

                    ///////////////////
                    // Tasklet code (T_l0_c0)
                    tmp_call_1_out = config_var_379_0_in_i_emiss_from_band_lw_0;
                    ///////////////////

                    tmp_call_1 = tmp_call_1_out;
                }

            }
        }

    }

    for (tmp_parfor_2_0 = 1; (tmp_parfor_2_0 <= __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0); tmp_parfor_2_0 = (tmp_parfor_2_0 + 1)) {
        {
            int tmp_index_38;
            int* v_config_var_379_i_band_from_reordered_g_lw;
            v_config_var_379_i_band_from_reordered_g_lw = (int*)(&(config_var_379->i_band_from_reordered_g_lw)[0]);

            {
                int config_var_379_0_in_i_band_from_reordered_g_lw_0 = v_config_var_379_i_band_from_reordered_g_lw[((- __f2dace_SOA_i_band_from_reordered_g_lw_d_0_s_8_config_0) + tmp_parfor_2_0)];
                int tmp_index_38_out;

                ///////////////////
                // Tasklet code (T_l1035_c1035)
                tmp_index_38_out = (config_var_379_0_in_i_band_from_reordered_g_lw_0 - __f2dace_SOA_i_emiss_from_band_lw_d_0_s_5_config_0);
                ///////////////////

                tmp_index_38 = tmp_index_38_out;
            }
            {
                int* config_var_379_0_in_i_emiss_from_band_lw = &v_config_var_379_i_emiss_from_band_lw[0];
                int tmp_index_38_0_in = tmp_index_38;
                double tmp_noncontig_2_out_0;

                ///////////////////
                // Tasklet code (T_l1035_c1035)
                tmp_noncontig_2_out_0 = config_var_379_0_in_i_emiss_from_band_lw[tmp_index_38_0_in];
                ///////////////////

                tmp_noncontig_2[(tmp_parfor_2_0 - 1)] = tmp_noncontig_2_out_0;
            }

        }

    }

    for (tmp_parfor_3_1 = 1; (tmp_parfor_3_1 <= __f2dace_SA_i_band_from_reordered_g_lw_d_0_s_8_config_0); tmp_parfor_3_1 = (tmp_parfor_3_1 + 1)) {
        for (tmp_parfor_2 = istartcol_var_380; (tmp_parfor_2 <= iendcol_var_381); tmp_parfor_2 = (tmp_parfor_2 + 1)) {

            tmp_index_39 = (tmp_parfor_2 - istartcol_var_380[0]);
            tmp_index_43 = (tmp_noncontig_2[(tmp_parfor_3_1 - 1)] - __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1);
            {
                double* v_this_var_378_lw_emissivity;
                v_this_var_378_lw_emissivity = (double*)(&(this_var_378->lw_emissivity)[0]);

                {
                    double this_var_378_0_in_lw_emissivity_0 = v_this_var_378_lw_emissivity[(((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * tmp_index_43) - __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_2)];
                    double tmp_noncontig_3_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    tmp_noncontig_3_out_0 = this_var_378_0_in_lw_emissivity_0;
                    ///////////////////

                    tmp_noncontig_3[(tmp_index_39 + ((tmp_parfor_3_1 - 1) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))] = tmp_noncontig_3_out_0;
                }

            }

        }

    }

    for (tmp_parfor_4 = __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1; (tmp_parfor_4 <= ((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 + __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) - 1)); tmp_parfor_4 = (tmp_parfor_4 + 1)) {
        for (tmp_parfor_3 = __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1; (tmp_parfor_3 <= ((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 + __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) - 1)); tmp_parfor_3 = (tmp_parfor_3 + 1)) {
            {

                {
                    double tmp_noncontig_3_0_in_0 = tmp_noncontig_3[(((- __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_3) + (((- __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) + tmp_parfor_4) * ((sym_iendcol_var_381 - sym_istartcol_var_380) + 1)))];
                    double tmp_libnode_0_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    tmp_libnode_0_out_0 = tmp_noncontig_3_0_in_0;
                    ///////////////////

                    tmp_libnode_0[(((__f2dace_SA_lw_emissivity_d_0_s_20_single_level_1 * ((- __f2dace_SOA_lw_emissivity_d_1_s_21_single_level_1) + tmp_parfor_4)) - __f2dace_SOA_lw_emissivity_d_0_s_20_single_level_1) + tmp_parfor_3)] = tmp_libnode_0_out_0;
                }

            }

        }

    }

    {

        transpose_sdfg_2_12_2(__state, &tmp_libnode_0[0], &tmp_call_3[0]);

    }

    for (tmp_parfor_6 = istartcol_var_380; (tmp_parfor_6 <= ((((sym_iendcol - sym_istartcol) + 1) + istartcol_var_380) - 1)); tmp_parfor_6 = (tmp_parfor_6 + 1)) {
        for (tmp_parfor_5 = 1; (tmp_parfor_5 <= 140); tmp_parfor_5 = (tmp_parfor_5 + 1)) {

            tmp_index_49 = (tmp_parfor_6 - istartcol_var_380[0]);
            tmp_index_51 = ((tmp_parfor_6 + (1 - istartcol_var_380[0])) - 1);
            {

                {
                    double tmp_call_3_0_in_0 = tmp_call_3[(((__f2dace_SA_lw_emissivity_d_1_s_21_single_level_1 * tmp_index_51) + tmp_parfor_5) - 1)];
                    double lw_albedo_var_382_out_0;

                    ///////////////////
                    // Tasklet code (T_l1035_c1035)
                    lw_albedo_var_382_out_0 = (1.0 - tmp_call_3_0_in_0);
                    ///////////////////

                    lw_albedo_var_382[(((140 * tmp_index_49) + tmp_parfor_5) - 1)] = lw_albedo_var_382_out_0;
                }

            }

        }

    }

    delete[] sw_albedo_band;
    delete[] tmp_call_3;
    delete[] tmp_libnode_0;
    delete[] tmp_noncontig_0;
    delete[] tmp_noncontig_1;
    delete[] tmp_noncontig_2;
    delete[] tmp_noncontig_3;
}

void __program_get_albedos_internal(get_albedos_state_t*__state, aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, double * __restrict__ lw_albedo_var_599, single_level_type* single_level, double * __restrict__ sw_albedo_diffuse_var_601, double * __restrict__ sw_albedo_direct_var_600, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    int nulerr;



    nulerr = 0;



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

        global_init_fn0_0_30_0(__state, nulout);

    }

    {

        get_albedos1_0_31_0(__state, config, &iendcol, &istartcol, single_level, &lw_albedo_var_599[0], &sw_albedo_diffuse_var_601[0], &sw_albedo_direct_var_600[0], sym_iendcol, iendcol, sym_istartcol, istartcol);

    }

}

DACE_EXPORTED void __program_get_albedos(get_albedos_state_t *__state, aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, double * __restrict__ lw_albedo_var_599, single_level_type* single_level, double * __restrict__ sw_albedo_diffuse_var_601, double * __restrict__ sw_albedo_direct_var_600, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    __program_get_albedos_internal(__state, aerosol, cloud, config, flux, gas, lw_albedo_var_599, single_level, sw_albedo_diffuse_var_601, sw_albedo_direct_var_600, thermodynamics, iendcol, istartcol, ncol, nlev, nulout, sym_iendcol, sym_istartcol);
}

DACE_EXPORTED get_albedos_state_t *__dace_init_get_albedos(aerosol_type* aerosol, cloud_type* cloud, config_type* config, flux_type* flux, gas_type* gas, double * __restrict__ lw_albedo_var_599, single_level_type* single_level, double * __restrict__ sw_albedo_diffuse_var_601, double * __restrict__ sw_albedo_direct_var_600, thermodynamics_type* thermodynamics, int iendcol, int istartcol, int ncol, int nlev, int nulout, int sym_iendcol, int sym_istartcol)
{
    int __result = 0;
    get_albedos_state_t *__state = new get_albedos_state_t;


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

DACE_EXPORTED int __dace_exit_get_albedos(get_albedos_state_t *__state)
{
    int __err = 0;
    delete __state;
    return __err;
}
