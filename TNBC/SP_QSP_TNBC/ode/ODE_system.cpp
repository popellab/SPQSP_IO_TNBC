#include "ODE_system.h"
    
#define SPVAR(x) NV_DATA_S(y)[x]
#define NSPVAR(x) ptrOde->_nonspecies_var[x]
#define PARAM(x) _class_parameter[x]
#define PFILE(x) param.getVal(x)

namespace CancerVCT{
#define QSP_W ODE_system::_QSP_weight

bool ODE_system::use_steady_state = false;
bool ODE_system::use_resection = false;
double ODE_system::_QSP_weight;  

ODE_system::ODE_system()
:CVODEBase()
{
    setupVariables();
    setupEvents();
    setupCVODE();
    update_y_other();
}

ODE_system::ODE_system(const ODE_system& c)
{
    setupCVODE();
}

ODE_system::~ODE_system()
{
}

void ODE_system::initSolver(realtype t){

    restore_y();
    int flag;

    flag = CVodeInit(_cvode_mem, f, t, _y);
    check_flag(&flag, "CVodeInit", 1);

    /* Call CVodeRootInit to specify the root function g */
    flag = CVodeRootInit(_cvode_mem, _nroot, g);
    check_flag(&flag, "CVodeRootInit", 1);
    
    	/*Do not do this. Event only trigger when turn from false to true.
	  If this is reset before trigger evaluation at the beginning of simulation,
	  t=0 events might be missed.*/
    //updateTriggerComponentConditionsOnValue(t);
    //resetEventTriggers();

    return;
}    

state_type ODE_system::_class_parameter = state_type(198, 0);

void ODE_system::setup_class_parameters(Param& param){
    //V_C, mw324aaeb5_4cad_4b2e_ae0c_8951b343415d, index: 0
    //Unit: metre^(3)
    _class_parameter[0] = PFILE(5) * 0.0010000000000000002;
    //V_P, mw1c783bd7_0eca_4d87_8829_58863541c89c, index: 1
    //Unit: metre^(3)
    _class_parameter[1] = PFILE(6) * 0.0010000000000000002;
    //V_LN, mwf7ed76cd_3f43_4158_88cb_5e0047c1285b, index: 2
    //Unit: metre^(3)
    _class_parameter[2] = PFILE(8) * 1.0000000000000013e-09;
    //V_e, mwfae11d7c_6828_4a44_b69a_64a0e3656d85, index: 3
    //Unit: metre^(3)
    _class_parameter[3] = PFILE(9) * 0.0010000000000000002;
    //A_e, mw8e40eb56_eaaa_4cf2_9b7b_ccf6f0f37bae, index: 4
    //Unit: metre^(2)
    _class_parameter[4] = PFILE(10) * 1e-12;
    //A_s, mwb87553d5_6aac_4f53_b8c4_a288e46766fa, index: 5
    //Unit: metre^(2)
    _class_parameter[5] = PFILE(11) * 1e-12;
    //syn_T_C1, mw15045a9f_a547_4e70_8547_035eaec325d3, index: 6
    //Unit: metre^(2)
    _class_parameter[6] = PFILE(12) * 1e-12;
    //syn_T_APC, mw8a3f7b63_2487_4017_931c_3f17eb64c5f6, index: 7
    //Unit: metre^(2)
    _class_parameter[7] = PFILE(13) * 1e-12;
    //k_cell_clear, mw7ded0f15_7079_4e08_8f41_f29d077c3fb3, index: 8
    //Unit: second^(-1)
    _class_parameter[8] = PFILE(134) * 1.15740740740741e-05;
    //cell, mw3219f4aa_c919_4aeb_a521_4e4e45145abb, index: 9
    //Unit: mole^(1)
    _class_parameter[9] = PFILE(135) * 1.66053872801495e-24;
    //day, mw873cf2a3_48bc_4d5b_b015_aec7bcaec40e, index: 10
    //Unit: second^(1)
    _class_parameter[10] = PFILE(136) * 86400.0;
    //vol_cell, mwf95c81c0_7393_4705_be71_106000d5f70a, index: 11
    //Unit: metre^(3)mole^(-1)
    _class_parameter[11] = PFILE(137) * 602214.1989999996;
    //vol_Tcell, mwfb1b43c2_b78e_4240_b335_b60a190f632a, index: 12
    //Unit: metre^(3)mole^(-1)
    _class_parameter[12] = PFILE(138) * 602214.1989999996;
    //V_Tmin, mwc3e516a2_b5ee_4123_b84b_53c342fe8800, index: 13
    //Unit: metre^(3)
    _class_parameter[13] = PFILE(139) * 1.0000000000000006e-06;
    //k_C1_growth, mw1c70db6c_706a_4c8d_ba5f_0f04df6321e4, index: 14
    //Unit: second^(-1)
    _class_parameter[14] = PFILE(153) * 1.15740740740741e-05;
    //C_max, mw42940e13_7ccc_462e_b8a7_56b54cf28b02, index: 15
    //Unit: mole^(1)
    _class_parameter[15] = PFILE(154) * 1.66053872801495e-24;
    //k_C1_death, mw172fb05a_325f_4eb3_989a_dddf14e6f7a6, index: 16
    //Unit: second^(-1)
    _class_parameter[16] = PFILE(155) * 1.15740740740741e-05;
    //k_C1_therapy, mw9ddd73b4_eaf9_41d5_9b15_d745457904d6, index: 17
    //Unit: second^(-1)
    _class_parameter[17] = PFILE(156) * 1.15740740740741e-05;
    //initial_tumour_diameter, mw09413063_6ae3_4dc8_bca8_6ea5b4c2c421, index: 18
    //Unit: metre^(1)
    _class_parameter[18] = PFILE(157) * 0.01;
    //n_T0_clones, mw40433098_1ce7_4592_a225_12c51791228e, index: 19
    //Unit: dimensionless^(1)
    _class_parameter[19] = PFILE(158) * 1.0;
    //Q_T0_in, mwe16f4e41_7017_43c6_adad_17fdf20c7e94, index: 20
    //Unit: mole^(1)second^(-1)
    _class_parameter[20] = PFILE(159) * 1.92191982409137e-29;
    //Q_T0_out, mwe97a7170_0f37_4057_8093_c11b2452cf24, index: 21
    //Unit: second^(-1)
    _class_parameter[21] = PFILE(160) * 1.15740740740741e-05;
    //k_T0_act, mw9a77ad14_cea3_4cf9_88e4_0f79f5a067b1, index: 22
    //Unit: second^(-1)
    _class_parameter[22] = PFILE(161) * 1.15740740740741e-05;
    //k_T0_pro, mw9688dcf2_ff84_489f_8d4e_e42e32628680, index: 23
    //Unit: second^(-1)
    _class_parameter[23] = PFILE(162) * 1.15740740740741e-05;
    //k_T0_death, mwc53d9795_0266_4b4a_9f2e_b1ee4a2215dd, index: 24
    //Unit: second^(-1)
    _class_parameter[24] = PFILE(163) * 1.15740740740741e-05;
    //q_T0_P_in, mw710ed142_a7e8_479d_af8a_06128d0736de, index: 25
    //Unit: second^(-1)
    _class_parameter[25] = PFILE(164) * 0.0166666666666667;
    //q_T0_P_out, mwfed80ebd_6b49_4181_ae55_82bf4a3821ef, index: 26
    //Unit: second^(-1)
    _class_parameter[26] = PFILE(165) * 1.15740740740741e-05;
    //q_T0_T_in, mwe71c2a74_e339_4309_b0a3_ffae4221d0f5, index: 27
    //Unit: metre^(-3)second^(-1)
    _class_parameter[27] = PFILE(166) * 16666.666666666693;
    //q_T0_LN_out, mwa6fb9cc8_60a7_4f9a_b0fe_54eaf246c14f, index: 28
    //Unit: second^(-1)
    _class_parameter[28] = PFILE(167) * 1.15740740740741e-05;
    //k_IL2_deg, mw20524179_3d46_4f4a_9273_f6c27e599e39, index: 29
    //Unit: second^(-1)
    _class_parameter[29] = PFILE(168) * 0.0166666666666667;
    //k_IL2_cons, mwcc08dc35_8b1e_47bd_a9f2_95cab1148c26, index: 30
    //Unit: second^(-1)
    _class_parameter[30] = PFILE(169) * 167281721944.444;
    //k_IL2_sec, mw8ce3daed_6e8f_48ae_a864_394ea94a5271, index: 31
    //Unit: second^(-1)
    _class_parameter[31] = PFILE(170) * 167281721944.444;
    //IL2_50, mw89746d66_728f_440c_9d81_fe1262546a6d, index: 32
    //Unit: metre^(-3)mole^(1)
    _class_parameter[32] = PFILE(171) * 1.0000000000000008e-06;
    //IL2_50_Treg, mwbfe0bc06_0327_4308_8218_4fc5c70e9eed, index: 33
    //Unit: metre^(-3)mole^(1)
    _class_parameter[33] = PFILE(172) * 1.0000000000000008e-06;
    //N0, mw47dba5cf_7969_4e24_9458_2f67cf42fe8a, index: 34
    //Unit: dimensionless^(1)
    _class_parameter[34] = PFILE(173) * 1.0;
    //N_costim, mw1fa7d989_eca1_4c20_a402_338a8e431153, index: 35
    //Unit: dimensionless^(1)
    _class_parameter[35] = PFILE(174) * 1.0;
    //N_IL2, mw2f8f7954_3483_42e8_853f_e3b8394d4e39, index: 36
    //Unit: dimensionless^(1)
    _class_parameter[36] = PFILE(175) * 1.0;
    //k_Treg, mw1f54ed76_4c5d_4e44_ab5d_06f95d6e6efa, index: 37
    //Unit: second^(-1)
    _class_parameter[37] = PFILE(176) * 1.15740740740741e-05;
    //n_T1_clones, mwc01ae599_f26d_4339_a279_dadb385cee85, index: 38
    //Unit: dimensionless^(1)
    _class_parameter[38] = PFILE(179) * 1.0;
    //Q_T1_in, mw44de9ada_c825_412a_9faf_07aee457ffad, index: 39
    //Unit: mole^(1)second^(-1)
    _class_parameter[39] = PFILE(180) * 1.92191982409137e-29;
    //Q_T1_out, mwbc50e857_6f5c_490b_b881_fc07604962cb, index: 40
    //Unit: second^(-1)
    _class_parameter[40] = PFILE(181) * 1.15740740740741e-05;
    //k_T1_act, mw85f39572_93f1_40d5_b76e_09a4cd2890a2, index: 41
    //Unit: second^(-1)
    _class_parameter[41] = PFILE(182) * 1.15740740740741e-05;
    //k_T1_pro, mwc735249e_b41f_42c3_aff5_bac0796441c7, index: 42
    //Unit: second^(-1)
    _class_parameter[42] = PFILE(183) * 1.15740740740741e-05;
    //k_T1_death, mwe0fd3e62_9009_4e5e_a5cd_92f7323fef48, index: 43
    //Unit: second^(-1)
    _class_parameter[43] = PFILE(184) * 1.15740740740741e-05;
    //q_T1_P_in, mwe257cda2_4621_475b_9ea2_b3fb17944ee1, index: 44
    //Unit: second^(-1)
    _class_parameter[44] = PFILE(185) * 0.0166666666666667;
    //q_T1_P_out, mw5b4f1d1c_4b8c_400a_98be_742b19edb605, index: 45
    //Unit: second^(-1)
    _class_parameter[45] = PFILE(186) * 1.15740740740741e-05;
    //q_T1_T_in, mw0998cc54_613a_4269_92f2_5ef45728a032, index: 46
    //Unit: metre^(-3)second^(-1)
    _class_parameter[46] = PFILE(187) * 16666.666666666693;
    //q_T1_LN_out, mw0f91adb1_e1a0_4b47_851e_57909aefce4e, index: 47
    //Unit: second^(-1)
    _class_parameter[47] = PFILE(188) * 1.15740740740741e-05;
    //k_T1, mwfad51e40_14e2_4b64_8721_2c6572fda68f, index: 48
    //Unit: second^(-1)
    _class_parameter[48] = PFILE(189) * 1.15740740740741e-05;
    //k_C_T1, mw584f9879_e233_41a1_9bb5_b5c96a1d5d64, index: 49
    //Unit: second^(-1)
    _class_parameter[49] = PFILE(190) * 1.15740740740741e-05;
    //k_APC_mat, mwb7b6477a_7e67_47bf_aa97_021b08e383da, index: 50
    //Unit: second^(-1)
    _class_parameter[50] = PFILE(192) * 1.15740740740741e-05;
    //k_APC_mig, mwa8981de5_f723_480d_b833_54b4f09e971f, index: 51
    //Unit: second^(-1)
    _class_parameter[51] = PFILE(193) * 1.15740740740741e-05;
    //k_APC_death, mw848c0b31_2892_4fd2_884f_93f3098a5cab, index: 52
    //Unit: second^(-1)
    _class_parameter[52] = PFILE(194) * 1.15740740740741e-05;
    //k_mAPC_death, mw43a9819d_5a51_4f4a_83fb_f7c7169d089c, index: 53
    //Unit: second^(-1)
    _class_parameter[53] = PFILE(195) * 1.15740740740741e-05;
    //APC0_T, mw76a9822b_e007_465e_8747_8986355dc5fc, index: 54
    //Unit: metre^(-3)mole^(1)
    _class_parameter[54] = PFILE(196) * 1.6605387280149534e-18;
    //APC0_LN, mw2d81ef83_faa2_4ab3_b98f_bbd25d901293, index: 55
    //Unit: metre^(-3)mole^(1)
    _class_parameter[55] = PFILE(197) * 1.6605387280149534e-18;
    //k_c, mw807e317e_560c_40f9_977c_794faaf32a55, index: 56
    //Unit: second^(-1)
    _class_parameter[56] = PFILE(198) * 1.15740740740741e-05;
    //c0, mw361830ff_0695_4180_9cba_01a7abac5a5c, index: 57
    //Unit: metre^(-3)mole^(1)
    _class_parameter[57] = PFILE(199) * 999.9999999999994;
    //c50, mw361a3d85_2712_4cb1_a85d_81e670724570, index: 58
    //Unit: metre^(-3)mole^(1)
    _class_parameter[58] = PFILE(200) * 999.9999999999994;
    //DAMPs, mw39344378_48ab_41b3_81c5_bd216d0ac262, index: 59
    //Unit: dimensionless^(1)
    _class_parameter[59] = PFILE(201) * 6.02214199e+23;
    //n_sites_APC, mw9b9ddae9_60d7_4217_b3a8_c593eed55a2e, index: 60
    //Unit: dimensionless^(1)
    _class_parameter[60] = PFILE(202) * 1.0;
    //kin, mw4334e045_9bd7_40f3_9892_a2bafa9975f1, index: 61
    //Unit: second^(-1)
    _class_parameter[61] = PFILE(203) * 1.15740740740741e-05;
    //kout, mw33340eab_dfa0_47dd_9f1c_cdc7e1b3f113, index: 62
    //Unit: second^(-1)
    _class_parameter[62] = PFILE(204) * 1.15740740740741e-05;
    //k_P0_up, mw3e4f357a_b51b_44e2_ad8a_2a0a244ce321, index: 63
    //Unit: mole^(-1)second^(-1)
    _class_parameter[63] = PFILE(205) * 6.97007174768519e+18;
    //k_xP0_deg, mw9861b3c5_d9b5_4765_a8f2_16e1bc5b7ceb, index: 64
    //Unit: second^(-1)
    _class_parameter[64] = PFILE(206) * 1.15740740740741e-05;
    //k_P0_deg, mwd38409eb_d3b3_4103_880a_a8a2f26cac12, index: 65
    //Unit: second^(-1)
    _class_parameter[65] = PFILE(207) * 1.15740740740741e-05;
    //k_p0_deg, mwf49d3257_f9f6_4ebd_89f3_50ed2ee044ab, index: 66
    //Unit: second^(-1)
    _class_parameter[66] = PFILE(208) * 1.15740740740741e-05;
    //k_P0_on, mwcedb38d5_362e_48c4_9105_4184cc176d83, index: 67
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[67] = PFILE(209) * 1.1574074074074112e-08;
    //k_P0_d1, mw8f784a6d_03b9_40e1_b462_f1ef3ef83525, index: 68
    //Unit: metre^(-3)mole^(1)
    _class_parameter[68] = PFILE(210) * 999.9999999999994;
    //p0_50, mwe0cecd41_c77b_4922_8d0f_3d5e52f5a07e, index: 69
    //Unit: metre^(-2)mole^(1)
    _class_parameter[69] = PFILE(211) * 1.66053872801495e-12;
    //P0_C1, mw40c735cd_cb8c_4a22_a177_d79985cd24a8, index: 70
    //Unit: metre^(-3)
    _class_parameter[70] = PFILE(212) * 6.022141989999979e+26;
    //A_syn, mw9d216aed_26e5_49e5_8658_05e01e5a50fe, index: 71
    //Unit: metre^(2)
    _class_parameter[71] = PFILE(213) * 1e-12;
    //A_Tcell, mw549b5ad2_5c81_413b_b9e1_60dfad846309, index: 72
    //Unit: metre^(2)
    _class_parameter[72] = PFILE(214) * 1e-12;
    //A_cell, mw56450c9f_9e06_49dd_aab7_30d28e024b5b, index: 73
    //Unit: metre^(2)
    _class_parameter[73] = PFILE(215) * 1e-12;
    //A_APC, mw76d68c82_d3fa_4139_9826_b768277bddd3, index: 74
    //Unit: metre^(2)
    _class_parameter[74] = PFILE(216) * 1e-12;
    //k_M1p0_TCR_on, mwb4548c6f_20fd_4bcc_9af2_ab9042a0b819, index: 75
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[75] = PFILE(217) * 602214199000.0;
    //k_M1p0_TCR_off, mw30b53315_debf_4bbb_a5d4_4069972ac311, index: 76
    //Unit: second^(-1)
    _class_parameter[76] = PFILE(218) * 1.0;
    //TCR_p0_tot, mw478e34f3_4cc6_42d8_87fa_426e1575d33b, index: 77
    //Unit: metre^(-2)mole^(1)
    _class_parameter[77] = PFILE(219) * 1.66053872801495e-12;
    //k_P1_up, mw0f37821a_f54f_4578_9448_d7bd24d5396a, index: 78
    //Unit: mole^(-1)second^(-1)
    _class_parameter[78] = PFILE(221) * 6.97007174768519e+18;
    //k_xP1_deg, mw488ac15f_bcda_4124_8550_362f2a3f11d1, index: 79
    //Unit: second^(-1)
    _class_parameter[79] = PFILE(222) * 1.15740740740741e-05;
    //k_P1_deg, mw390cedeb_b602_4a82_b65a_77ee6c205e65, index: 80
    //Unit: second^(-1)
    _class_parameter[80] = PFILE(223) * 1.15740740740741e-05;
    //k_p1_deg, mw7b6465a9_d9ae_41db_b2c8_90898bba0aa1, index: 81
    //Unit: second^(-1)
    _class_parameter[81] = PFILE(224) * 1.15740740740741e-05;
    //k_P1_on, mwec31aa0c_e8c5_4de8_831a_33d8f515277d, index: 82
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[82] = PFILE(225) * 1.1574074074074112e-08;
    //k_P1_d1, mw096496aa_995d_4e70_8ccf_53cce7c1ce63, index: 83
    //Unit: metre^(-3)mole^(1)
    _class_parameter[83] = PFILE(226) * 999.9999999999994;
    //p1_50, mwe8173793_2ca6_4a02_a6a2_c18ef3a1c6ed, index: 84
    //Unit: metre^(-2)mole^(1)
    _class_parameter[84] = PFILE(227) * 1.66053872801495e-12;
    //P1_C1, mw8651a4a0_5619_494d_ba98_151c9760efea, index: 85
    //Unit: metre^(-3)
    _class_parameter[85] = PFILE(228) * 6.022141989999979e+26;
    //k_M1p1_TCR_on, mw9d2808da_3bba_4c89_92e3_43c55d47c29c, index: 86
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[86] = PFILE(229) * 602214199000.0;
    //k_M1p1_TCR_off, mwa302dee9_68d9_45e3_9804_fd3ce984d504, index: 87
    //Unit: second^(-1)
    _class_parameter[87] = PFILE(230) * 1.0;
    //k_M1p1_TCR_p, mw0aecb575_312d_44d8_8725_e3b6db206757, index: 88
    //Unit: second^(-1)
    _class_parameter[88] = PFILE(231) * 1.0;
    //phi_M1p1_TCR, mw99631140_328b_4a2d_b58c_28217a130b59, index: 89
    //Unit: second^(-1)
    _class_parameter[89] = PFILE(232) * 1.0;
    //N_M1p1_TCR, mwac8e27a0_ada1_448c_8e9b_d68ff3042555, index: 90
    //Unit: dimensionless^(1)
    _class_parameter[90] = PFILE(233) * 1.0;
    //TCR_p1_tot, mwb9ca1eb2_bf52_407a_99e2_ad0db4ea46d2, index: 91
    //Unit: metre^(-2)mole^(1)
    _class_parameter[91] = PFILE(234) * 1.66053872801495e-12;
    //kon_PD1_PDL1, mw4afa9502_bda0_4053_ab51_295aeaf6d4e5, index: 92
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[92] = PFILE(236) * 1000000000000.0;
    //q_P_nivo, mw60e3a0dd_960e_4b0f_b930_bc5971f1f1d6, index: 93
    //Unit: metre^(3)second^(-1)
    _class_parameter[93] = PFILE(237) * 0.0010000000000000007;
    //q_T_nivo, mw5a7da654_87e0_42bd_962d_3f398e18eae3, index: 94
    //Unit: metre^(3)second^(-1)
    _class_parameter[94] = PFILE(238) * 1.0000000000000006e-06;
    //q_LN_nivo, mwb7096f53_5701_44cd_b45b_590a8bb1ea0e, index: 95
    //Unit: metre^(3)second^(-1)
    _class_parameter[95] = PFILE(239) * 1.0000000000000006e-06;
    //q_LD_nivo, mw545d39cf_90c4_4fdb_852a_881f4b415bc8, index: 96
    //Unit: second^(-1)
    _class_parameter[96] = PFILE(240) * 0.0166666666666667;
    //k_cl_nivo, mw4bf315b4_4448_4780_8584_1b7140224377, index: 97
    //Unit: metre^(3)second^(-1)
    _class_parameter[97] = PFILE(241) * 1.1574074074074112e-08;
    //gamma_C_nivo, mw21275fa2_2557_4c32_a3d9_e665cba46851, index: 98
    //Unit: dimensionless^(1)
    _class_parameter[98] = PFILE(242) * 1.0;
    //gamma_P_nivo, mw2aa4543f_6511_4d5a_8ae0_37941c8ac97a, index: 99
    //Unit: dimensionless^(1)
    _class_parameter[99] = PFILE(243) * 1.0;
    //gamma_T_nivo, mw045b51d2_15ee_4291_98b2_b840b11bdbb1, index: 100
    //Unit: dimensionless^(1)
    _class_parameter[100] = PFILE(244) * 1.0;
    //gamma_LN_nivo, mw22861acd_d08d_4998_9ba0_66bebad3d43d, index: 101
    //Unit: dimensionless^(1)
    _class_parameter[101] = PFILE(245) * 1.0;
    //q_P_durv, mw01262d82_a03b_4ff1_ba2d_140a2f9b9077, index: 102
    //Unit: metre^(3)second^(-1)
    _class_parameter[102] = PFILE(246) * 0.0010000000000000007;
    //q_T_durv, mw8e9da096_c128_4255_af52_242692465df2, index: 103
    //Unit: metre^(3)second^(-1)
    _class_parameter[103] = PFILE(247) * 1.0000000000000006e-06;
    //q_LN_durv, mw22d415c5_0549_4274_9dc3_2c1384032c83, index: 104
    //Unit: metre^(3)second^(-1)
    _class_parameter[104] = PFILE(248) * 1.0000000000000006e-06;
    //q_LD_durv, mw444563e6_c678_411f_9f59_b954b8ac90d1, index: 105
    //Unit: second^(-1)
    _class_parameter[105] = PFILE(249) * 0.0166666666666667;
    //k_cl_durv, mw5f492e39_23ce_436b_802b_e00c5e1af9aa, index: 106
    //Unit: metre^(3)second^(-1)
    _class_parameter[106] = PFILE(250) * 1.1574074074074112e-08;
    //gamma_C_durv, mw7c74c8ea_18fa_452f_b63b_761dbd344c95, index: 107
    //Unit: dimensionless^(1)
    _class_parameter[107] = PFILE(251) * 1.0;
    //gamma_P_durv, mwf37130b8_2bd1_4860_b285_7b7ec2c8cfef, index: 108
    //Unit: dimensionless^(1)
    _class_parameter[108] = PFILE(252) * 1.0;
    //gamma_T_durv, mw78d61efb_1b92_4aef_9987_bcc0a69e45d8, index: 109
    //Unit: dimensionless^(1)
    _class_parameter[109] = PFILE(253) * 1.0;
    //gamma_LN_durv, mw251c5b11_25b6_4a87_82f4_540412af5ea7, index: 110
    //Unit: dimensionless^(1)
    _class_parameter[110] = PFILE(254) * 1.0;
    //q_P_ipi, mw022c8e8e_5eb7_4257_90d2_3b67c1acc78a, index: 111
    //Unit: metre^(3)second^(-1)
    _class_parameter[111] = PFILE(255) * 0.0010000000000000007;
    //q_T_ipi, mw51f039e0_edd8_4078_847d_bd1d134f26b3, index: 112
    //Unit: metre^(3)second^(-1)
    _class_parameter[112] = PFILE(256) * 1.0000000000000006e-06;
    //q_LN_ipi, mw361e65bd_36df_45a9_9d31_797b502e79df, index: 113
    //Unit: metre^(3)second^(-1)
    _class_parameter[113] = PFILE(257) * 1.0000000000000006e-06;
    //q_LD_ipi, mw81f0ae42_35dc_4d37_9511_6440adf8785e, index: 114
    //Unit: second^(-1)
    _class_parameter[114] = PFILE(258) * 0.0166666666666667;
    //k_cl_ipi, mwac77271e_b433_49c9_a2b4_66f9c6afe0a2, index: 115
    //Unit: metre^(3)second^(-1)
    _class_parameter[115] = PFILE(259) * 1.1574074074074112e-08;
    //gamma_C_ipi, mwfd4bac79_e449_4821_98e4_88759688637e, index: 116
    //Unit: dimensionless^(1)
    _class_parameter[116] = PFILE(260) * 1.0;
    //gamma_P_ipi, mw85fb5088_d435_4361_8dd2_f157750e5144, index: 117
    //Unit: dimensionless^(1)
    _class_parameter[117] = PFILE(261) * 1.0;
    //gamma_T_ipi, mw1490e1c2_26d1_4b83_9cf3_16210cf274a9, index: 118
    //Unit: dimensionless^(1)
    _class_parameter[118] = PFILE(262) * 1.0;
    //gamma_LN_ipi, mwa30cd3bf_494b_4af8_9fc9_6c8740e45e1e, index: 119
    //Unit: dimensionless^(1)
    _class_parameter[119] = PFILE(263) * 1.0;
    //kon_PD1_PDL2, mw2810ffba_871f_41e3_b13f_3770a16a9eb7, index: 120
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[120] = PFILE(264) * 1000000000000.0;
    //kon_PD1_nivo, mw52cfa457_aa56_45dd_9f67_234d3ca5fa8f, index: 121
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[121] = PFILE(265) * 0.0010000000000000007;
    //kon_PDL1_durv, mwfd3f3a6d_09e0_46f6_a1ef_86b8ed405d4a, index: 122
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[122] = PFILE(266) * 0.0010000000000000007;
    //kon_CD28_CD80, mwcf7cc6fd_070d_4c39_b57b_d9562a168ad6, index: 123
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[123] = PFILE(267) * 1000000000000.0;
    //kon_CD28_CD86, mw2e05d635_5432_4b41_b767_7c76be2f9bdc, index: 124
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[124] = PFILE(268) * 1000000000000.0;
    //kon_CTLA4_CD80, mw31f48165_cf74_4d96_9bd3_211f9c65a8d5, index: 125
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[125] = PFILE(269) * 1000000000000.0;
    //kon_CTLA4_CD86, mwebab6bc5_f1e2_45f5_83aa_0d7071343f18, index: 126
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[126] = PFILE(270) * 1000000000000.0;
    //kon_CD80_PDL1, mw51d56879_f327_41c4_a98b_5d8197b725eb, index: 127
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[127] = PFILE(271) * 1000000000000.0;
    //kon_CTLA4_ipi, mw9c833d89_dc82_4a2d_adde_e58b86191d37, index: 128
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[128] = PFILE(272) * 0.0010000000000000007;
    //koff_PD1_PDL1, mw9c199711_6bf5_4d86_b86f_840eb4d89b8b, index: 129
    //Unit: second^(-1)
    _class_parameter[129] = PFILE(273) * 1.0;
    //koff_PD1_PDL2, mwee8ce984_4a3d_4402_876f_db1e59dfc831, index: 130
    //Unit: second^(-1)
    _class_parameter[130] = PFILE(274) * 1.0;
    //koff_PD1_nivo, mwa115f6fc_b5ae_4a7e_94bb_14dee2f59556, index: 131
    //Unit: second^(-1)
    _class_parameter[131] = PFILE(275) * 1.0;
    //koff_PDL1_durv, mw6e44bf8a_4306_4beb_b1ff_25d04564dbfa, index: 132
    //Unit: second^(-1)
    _class_parameter[132] = PFILE(276) * 1.0;
    //koff_CD28_CD80, mw621901cc_f50e_4bc7_adf8_b5c351fb09a8, index: 133
    //Unit: second^(-1)
    _class_parameter[133] = PFILE(277) * 1.0;
    //koff_CD28_CD86, mw28552336_0c85_437c_9998_2e8f9ce68de9, index: 134
    //Unit: second^(-1)
    _class_parameter[134] = PFILE(278) * 1.0;
    //koff_CTLA4_CD80, mw84e2deb9_1cd3_450c_b2b8_94c43d2162d7, index: 135
    //Unit: second^(-1)
    _class_parameter[135] = PFILE(279) * 1.0;
    //koff_CTLA4_CD86, mw6ac546f0_4077_45ce_9e73_24328166614d, index: 136
    //Unit: second^(-1)
    _class_parameter[136] = PFILE(280) * 1.0;
    //koff_CD80_PDL1, mwd7ddebf0_8a9e_4baa_86c3_baac91692dbd, index: 137
    //Unit: second^(-1)
    _class_parameter[137] = PFILE(281) * 1.0;
    //koff_CTLA4_ipi, mwae815775_3bd1_43bc_bda3_5385097b2101, index: 138
    //Unit: second^(-1)
    _class_parameter[138] = PFILE(282) * 1.0;
    //Chi_PD1_nivo, mw0281e875_d54f_415a_abed_8a696977a6ac, index: 139
    //Unit: metre^(-1)
    _class_parameter[139] = PFILE(283) * 999999999.9999999;
    //Chi_PDL1_durv, mw23919693_9302_4284_811c_aab220316fcf, index: 140
    //Unit: metre^(-1)
    _class_parameter[140] = PFILE(284) * 999999999.9999999;
    //Chi_CTLA4_ipi, mw219f39b9_6590_4919_b37b_ded9dc041530, index: 141
    //Unit: metre^(-1)
    _class_parameter[141] = PFILE(285) * 999999999.9999999;
    //PD1_50, mw127481e6_dc43_448d_82a3_46d42b7066eb, index: 142
    //Unit: metre^(-2)mole^(1)
    _class_parameter[142] = PFILE(286) * 1.66053872801495e-12;
    //n_PD1, mwd9cae7f6_a1ac_4796_8932_3f5f5a48c324, index: 143
    //Unit: dimensionless^(1)
    _class_parameter[143] = PFILE(287) * 1.0;
    //CD28_CD8X_50, mwced7b380_945f_4154_9aca_b26abb9cb272, index: 144
    //Unit: metre^(-2)mole^(1)
    _class_parameter[144] = PFILE(288) * 1.66053872801495e-12;
    //n_CD28_CD8X, mwe53b5a8c_224b_4bcf_81ba_c47c32c41330, index: 145
    //Unit: dimensionless^(1)
    _class_parameter[145] = PFILE(289) * 1.0;
    //T_PD1_total, mw20938ee3_6910_4dcd_9b63_65a7ea804154, index: 146
    //Unit: mole^(1)
    _class_parameter[146] = PFILE(290) * 1.66053872801495e-24;
    //T_CD28_total, mw49777903_ed7f_4724_874a_ed3187284328, index: 147
    //Unit: mole^(1)
    _class_parameter[147] = PFILE(291) * 1.66053872801495e-24;
    //T_CTLA4_syn, mw20a050c8_619e_4dfe_84c3_29ca311df2bd, index: 148
    //Unit: mole^(1)
    _class_parameter[148] = PFILE(292) * 1.66053872801495e-24;
    //T_PDL1_total, mwe558348e_bb4a_41eb_8d7c_4b9357fe9e0b, index: 149
    //Unit: mole^(1)
    _class_parameter[149] = PFILE(293) * 1.66053872801495e-24;
    //C1_PDL1_total, mw8bd19d82_a32e_45ef_9d8f_6aa5ca8a93f9, index: 150
    //Unit: mole^(1)
    _class_parameter[150] = PFILE(294) * 1.66053872801495e-24;
    //C1_PDL2_total, mw38f890ad_f9dd_462a_a955_5d9acea9319b, index: 151
    //Unit: mole^(1)
    _class_parameter[151] = PFILE(295) * 1.66053872801495e-24;
    //C1_CD80_total, mw39923c4b_3d00_4cb4_9a3b_ef54f508bae4, index: 152
    //Unit: mole^(1)
    _class_parameter[152] = PFILE(296) * 1.66053872801495e-24;
    //C1_CD86_total, mw23acea19_3932_4b79_8f99_c8f34fc7bfae, index: 153
    //Unit: mole^(1)
    _class_parameter[153] = PFILE(297) * 1.66053872801495e-24;
    //APC_PDL1_total, mw4bd623f6_4a25_48d3_af3e_43a663a650e1, index: 154
    //Unit: mole^(1)
    _class_parameter[154] = PFILE(298) * 1.66053872801495e-24;
    //APC_PDL2_total, mw3e8d5b0e_8314_412b_a4c9_284487ad784b, index: 155
    //Unit: mole^(1)
    _class_parameter[155] = PFILE(299) * 1.66053872801495e-24;
    //APC_CD80_total, mwb5ae97a4_eff7_4514_9e35_21ae08463370, index: 156
    //Unit: mole^(1)
    _class_parameter[156] = PFILE(300) * 1.66053872801495e-24;
    //APC_CD86_total, mw6b465bc7_064f_42b9_89e2_6599d1323631, index: 157
    //Unit: mole^(1)
    _class_parameter[157] = PFILE(301) * 1.66053872801495e-24;
    //Treg_CTLA4_tot, mw17d4c2e3_8567_4775_bf1a_8a4fef2bf566, index: 158
    //Unit: mole^(1)
    _class_parameter[158] = PFILE(302) * 1.66053872801495e-24;
    //Treg_CTLA4_50, mw58aefda1_8d32_4217_b551_be1c2dcf01f2, index: 159
    //Unit: mole^(1)
    _class_parameter[159] = PFILE(303) * 1.66053872801495e-24;
    //n_Treg_CTLA4, mw10b74df2_4762_42df_907c_3626b5b9a9cf, index: 160
    //Unit: dimensionless^(1)
    _class_parameter[160] = PFILE(304) * 1.0;
    //k_CTLA4_ADCC, mw98ee6307_028e_47f7_a57a_216dbc086479, index: 161
    //Unit: second^(-1)
    _class_parameter[161] = PFILE(305) * 1.15740740740741e-05;
    //k_rec_MDSC, mw1c763d18_7b18_458f_9683_1a90956d47a8, index: 162
    //Unit: second^(-1)
    _class_parameter[162] = PFILE(308) * 1.15740740740741e-05;
    //kd_MDSC, mw59db0e06_6747_4f7e_89bc_c6a5cb4a59c0, index: 163
    //Unit: second^(-1)
    _class_parameter[163] = PFILE(309) * 1.15740740740741e-05;
    //IC50_ENT_C, mwb5579903_cf2c_463b_8ea3_505bf50feab9, index: 164
    //Unit: metre^(-3)mole^(1)
    _class_parameter[164] = PFILE(310) * 999.9999999999994;
    //k_deg_CCL2, mw77c59a82_1821_48fa_b3d9_b175b1bb8a84, index: 165
    //Unit: second^(-1)
    _class_parameter[165] = PFILE(311) * 0.000277777777777778;
    //k_deg_NO, mw79c59ffb_f812_491f_af32_4cec33b7a4aa, index: 166
    //Unit: second^(-1)
    _class_parameter[166] = PFILE(312) * 1.15740740740741e-05;
    //k_deg_ArgI, mw91a23e45_f98f_411c_8e55_420464499b3a, index: 167
    //Unit: second^(-1)
    _class_parameter[167] = PFILE(313) * 1.15740740740741e-05;
    //k_sec_CCL2, mw38482bc8_1030_4e22_8d96_d5e1c98aab37, index: 168
    //Unit: second^(-1)
    _class_parameter[168] = PFILE(314) * 6970071747.68519;
    //k_sec_NO, mw768f4bd3_93b4_4c9a_a00f_61aa7d73a3f0, index: 169
    //Unit: second^(-1)
    _class_parameter[169] = PFILE(315) * 6970071747.68519;
    //k_sec_ArgI, mwb074dca8_ed98_4725_8ddf_e8b63c9ef57b, index: 170
    //Unit: second^(-1)
    _class_parameter[170] = PFILE(316) * 6970071747685.19;
    //IC50_ENT_NO, mw6a119ae9_218d_4d8f_a16c_db3a6b36ab1f, index: 171
    //Unit: metre^(-3)mole^(1)
    _class_parameter[171] = PFILE(317) * 999.9999999999994;
    //ki_Treg, mwb3c44fa0_a39d_4ac6_9ead_f3c42c74ab43, index: 172
    //Unit: second^(-1)
    _class_parameter[172] = PFILE(318) * 1.15740740740741e-05;
    //IC50_ArgI_CTL, mwd3bb2d6d_c3da_431a_aa2a_b5c514c805c2, index: 173
    //Unit: metre^(-3)mole^(1)
    _class_parameter[173] = PFILE(319) * 999.9999999999994;
    //IC50_NO_CTL, mw893460ad_0455_463e_b617_7ab39b0a50e6, index: 174
    //Unit: metre^(-3)mole^(1)
    _class_parameter[174] = PFILE(320) * 999.9999999999994;
    //EC50_CCL2_rec, mw499fe725_055e_4b57_8ec0_792073e47471, index: 175
    //Unit: metre^(-3)mole^(1)
    _class_parameter[175] = PFILE(321) * 999.9999999999994;
    //EC50_ArgI_Treg, mwfcf2c43c_060e_43c7_8b0b_44bd977df1c7, index: 176
    //Unit: metre^(-3)mole^(1)
    _class_parameter[176] = PFILE(322) * 999.9999999999994;
    //MDSC_max, mw2ecfcb93_fd1c_4d00_85f3_72b6af6c4c14, index: 177
    //Unit: metre^(-3)mole^(1)
    _class_parameter[177] = PFILE(323) * 1.6605387280149534e-18;
    //Treg_max, mw9bccc5cf_c721_46fa_8585_dbe934728213, index: 178
    //Unit: metre^(-3)mole^(1)
    _class_parameter[178] = PFILE(324) * 1.6605387280149534e-18;
    //IC50_ENT_CCL2, mw23257c95_bdbd_4464_8599_a174ca490f25, index: 179
    //Unit: metre^(-3)mole^(1)
    _class_parameter[179] = PFILE(325) * 999.9999999999994;
    //k_brec_MDSC, mw7520abb5_5e88_4318_9002_236b1c23c5a4, index: 180
    //Unit: second^(-1)
    _class_parameter[180] = PFILE(326) * 1.15740740740741e-05;
    //IC50_ENT_ArgI, mw92419197_a402_4a6a_8ab2_52c8b8227c41, index: 181
    //Unit: metre^(-3)mole^(1)
    _class_parameter[181] = PFILE(327) * 999.9999999999994;
    //k_a1_ENT, mw645820bf_d861_4b50_b748_41a071cf3be8, index: 182
    //Unit: second^(-1)
    _class_parameter[182] = PFILE(328) * 0.000277777777777778;
    //k_a2_ENT, mw297d24d7_fd90_4655_9f0a_4715e7af219c, index: 183
    //Unit: second^(-1)
    _class_parameter[183] = PFILE(329) * 0.000277777777777778;
    //k_cln_ENT, mw6584e1a9_1a7c_4e50_8469_f9eef8351812, index: 184
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[184] = PFILE(330) * 0.277777777777778;
    //Kc_ENT, mw261f364e_ef49_473b_ae7b_573069373b4d, index: 185
    //Unit: metre^(-3)mole^(1)
    _class_parameter[185] = PFILE(331) * 999.9999999999994;
    //lagP, mw99549309_3c3e_483b_8281_987153ea6944, index: 186
    //Unit: second^(1)
    _class_parameter[186] = PFILE(332) * 3600.0;
    //durP, mwc3ac59ca_c167_489c_87ea_f08839822c70, index: 187
    //Unit: second^(1)
    _class_parameter[187] = PFILE(333) * 3600.0;
    //k_dose2, mw0fdc8f2e_d35f_47f2_b152_945656feb63c, index: 188
    //Unit: second^(-1)
    _class_parameter[188] = PFILE(334) * 0.000277777777777778;
    //q_P_ENT, mw5600406b_3913_4120_82bc_3f6b9faf450d, index: 189
    //Unit: metre^(3)second^(-1)
    _class_parameter[189] = PFILE(335) * 1.0000000000000006e-06;
    //q_T_ENT, mw3ccc3fde_844f_4e71_9c1d_572eaad92b85, index: 190
    //Unit: metre^(3)second^(-1)
    _class_parameter[190] = PFILE(336) * 1.0000000000000006e-06;
    //q_LN_ENT, mw5e741e02_4403_4d5b_be38_ad84476ad587, index: 191
    //Unit: metre^(3)second^(-1)
    _class_parameter[191] = PFILE(337) * 1.0000000000000006e-06;
    //q_LD_ENT, mw0ffd106e_a126_4b2c_9854_33c61ca94d48, index: 192
    //Unit: second^(-1)
    _class_parameter[192] = PFILE(338) * 0.0166666666666667;
    //k_cl_ENT, mw7d1c006e_bdc2_4e67_829a_d9d7ba4e3744, index: 193
    //Unit: second^(-1)
    _class_parameter[193] = PFILE(339) * 0.000277777777777778;
    //gamma_C_ENT, mw6a630837_928f_4afc_be3c_badd75662da9, index: 194
    //Unit: dimensionless^(1)
    _class_parameter[194] = PFILE(340) * 1.0;
    //gamma_P_ENT, mw9895ba75_2182_4158_9dea_94c6f63f613d, index: 195
    //Unit: dimensionless^(1)
    _class_parameter[195] = PFILE(341) * 1.0;
    //gamma_T_ENT, mw02abe3a0_713f_487e_9a0c_fee463c4ec0c, index: 196
    //Unit: dimensionless^(1)
    _class_parameter[196] = PFILE(342) * 1.0;
    //gamma_LN_ENT, mw330df805_5783_43c0_9ad1_133b1655014d, index: 197
    //Unit: dimensionless^(1)
    _class_parameter[197] = PFILE(343) * 1.0;
}

void ODE_system::setupVariables(void){

    _species_var = std::vector<realtype>(120, 0);
    _nonspecies_var = std::vector<realtype>(0, 0);
    //species not part of ode left-hand side
    _species_other =  std::vector<realtype>(0, 0);
    
    return;
}


void ODE_system::setup_instance_varaibles(Param& param){

    //V_C.T0, mw69f7efe7_ac08_48bd_a951_50b4fbfab6eb, index: 0
    //Unit: mole^(1)
    _species_var[0] = PFILE(14) * 1.66053872801495e-24;
    //V_C.T1, mw8b7e1ba3_459e_43bf_acd4_ca0b2b71a72c, index: 1
    //Unit: mole^(1)
    _species_var[1] = PFILE(15) * 1.66053872801495e-24;
    //V_C.nivo, mw092b1fb2_83fd_4769_ac98_cf59f0bb2784, index: 2
    //Unit: mole^(1)metre^(-3)
    _species_var[2] = PFILE(16) * 999.9999999999999;
    //V_C.durv, mwf124d3c7_c221_45c2_9a79_cead823b24ee, index: 3
    //Unit: mole^(1)metre^(-3)
    _species_var[3] = PFILE(17) * 999.9999999999999;
    //V_C.ipi, mw444e73d2_30fb_4f3f_836e_05322210264a, index: 4
    //Unit: mole^(1)metre^(-3)
    _species_var[4] = PFILE(18) * 999.9999999999999;
    //V_C.ENT, mw0a530416_36c0_4c4c_99db_b0e03f18805b, index: 5
    //Unit: mole^(1)metre^(-3)
    _species_var[5] = PFILE(19) * 999.9999999999999;
    //V_C.ENT_Buccal, mw115206c3_1b63_42f2_a4c3_0f04995e2a65, index: 6
    //Unit: mole^(1)metre^(-3)
    _species_var[6] = PFILE(20) * 999.9999999999999;
    //V_C.ENT_GI, mw4c0976aa_d375_4932_9573_ed585072c3dc, index: 7
    //Unit: mole^(1)metre^(-3)
    _species_var[7] = PFILE(21) * 999.9999999999999;
    //V_C.Dose2, mwd6e1b481_4751_4d04_ac02_5fb6a8c321ec, index: 8
    //Unit: mole^(1)metre^(-3)
    _species_var[8] = PFILE(22) * 999.9999999999999;
    //V_P.T0, mwb8d31ba6_f2b3_422b_826d_ed90bd4790f2, index: 9
    //Unit: mole^(1)
    _species_var[9] = PFILE(23) * 1.66053872801495e-24;
    //V_P.T1, mwb4354006_cde2_457e_8b5a_8019a03cfa19, index: 10
    //Unit: mole^(1)
    _species_var[10] = PFILE(24) * 1.66053872801495e-24;
    //V_P.nivo, mw5fe8db7e_a2bc_4869_afa7_7f5ca0145daf, index: 11
    //Unit: mole^(1)metre^(-3)
    _species_var[11] = PFILE(25) * 999.9999999999999;
    //V_P.durv, mw509329f2_a9f6_4506_baaf_1511256fcb86, index: 12
    //Unit: mole^(1)metre^(-3)
    _species_var[12] = PFILE(26) * 999.9999999999999;
    //V_P.ipi, mw664b5eef_723d_49ec_b1d9_acc99bd3d148, index: 13
    //Unit: mole^(1)metre^(-3)
    _species_var[13] = PFILE(27) * 999.9999999999999;
    //V_P.Treg_CTLA4, mwd3fce12a_1681_487d_bd43_3f36da38db3a, index: 14
    //Unit: mole^(1)
    _species_var[14] = PFILE(28) * 1.66053872801495e-24;
    //V_P.Treg_CTLA4_ipi, mw3415a9cc_b566_4856_aaef_4d86f91f703b, index: 15
    //Unit: mole^(1)
    _species_var[15] = PFILE(29) * 1.66053872801495e-24;
    //V_P.Treg_CTLA4_ipi_CTLA4, mw29b7a526_7019_4266_9931_9827a36a92a9, index: 16
    //Unit: mole^(1)
    _species_var[16] = PFILE(30) * 1.66053872801495e-24;
    //V_P.ENT, mwdab736c4_6bbd_4ddb_9748_6fbc4abc0c23, index: 17
    //Unit: mole^(1)metre^(-3)
    _species_var[17] = PFILE(31) * 999.9999999999999;
    //V_T.C_x, mw3fb8ca45_3a74_4e8a_8db4_87f50476daca, index: 18
    //Unit: mole^(1)
    _species_var[18] = PFILE(32) * 1.66053872801495e-24;
    //V_T.T_exh, mwce590c3a_6cf5_45ba_8bb2_904a2e3dde91, index: 19
    //Unit: mole^(1)
    _species_var[19] = PFILE(33) * 1.66053872801495e-24;
    //V_T.C1, mwe7ebcac4_fadd_4936_8667_c34cf469fe3a, index: 20
    //Unit: mole^(1)
    _species_var[20] = PFILE(34) * 1.66053872801495e-24;
    //V_T.T0, mwddfbf0aa_a5c3_4212_8a86_85208843c43f, index: 21
    //Unit: mole^(1)
    _species_var[21] = PFILE(35) * 1.66053872801495e-24;
    //V_T.T1, mw29266435_9932_4c05_9b63_03d1754056f1, index: 22
    //Unit: mole^(1)
    _species_var[22] = PFILE(36) * 1.66053872801495e-24;
    //V_T.APC, mw73386c66_d30c_401e_98d0_eccb28102e1a, index: 23
    //Unit: mole^(1)
    _species_var[23] = PFILE(37) * 1.66053872801495e-24;
    //V_T.mAPC, mw58612907_55a0_4a82_8cf2_f75a5c495c0b, index: 24
    //Unit: mole^(1)
    _species_var[24] = PFILE(38) * 1.66053872801495e-24;
    //V_T.c, mw6ac93b40_ddaa_4388_b5e6_f541274b93df, index: 25
    //Unit: mole^(1)metre^(-3)
    _species_var[25] = PFILE(39) * 1000000.0;
    //V_T.nivo, mw4e5b59d9_0b34_42a8_95c4_29a1da1a3161, index: 26
    //Unit: mole^(1)metre^(-3)
    _species_var[26] = PFILE(40) * 1000000.0;
    //V_T.durv, mwfc65be48_9503_4263_9409_c76d8526daa1, index: 27
    //Unit: mole^(1)metre^(-3)
    _species_var[27] = PFILE(41) * 1000000.0;
    //V_T.ipi, mwdab20b92_e602_424b_9dbc_1cc52e0946d6, index: 28
    //Unit: mole^(1)metre^(-3)
    _species_var[28] = PFILE(42) * 1000000.0;
    //V_T.Treg_CTLA4, mweb4943a4_297b_4316_8595_07d7fe602cb8, index: 29
    //Unit: mole^(1)
    _species_var[29] = PFILE(43) * 1.66053872801495e-24;
    //V_T.Treg_CTLA4_ipi, mwea43ca96_26b6_408a_8d84_0539f63fc11f, index: 30
    //Unit: mole^(1)
    _species_var[30] = PFILE(44) * 1.66053872801495e-24;
    //V_T.Treg_CTLA4_ipi_CTLA4, mwf4a7b54b_9a68_48f3_91c2_6d25453adb85, index: 31
    //Unit: mole^(1)
    _species_var[31] = PFILE(45) * 1.66053872801495e-24;
    //V_T.MDSC, mw83cbb11d_d1ab_4c5d_ab9b_31621cc6e61a, index: 32
    //Unit: mole^(1)
    _species_var[32] = PFILE(46) * 1.66053872801495e-24;
    //V_T.CCL2, mwabb84d39_4d5e_43df_b89a_95c16ee105f4, index: 33
    //Unit: mole^(1)metre^(-3)
    _species_var[33] = PFILE(47) * 1000000.0;
    //V_T.NO, mwa8265122_d899_4576_a773_f4e09670d29c, index: 34
    //Unit: mole^(1)metre^(-3)
    _species_var[34] = PFILE(48) * 1000000.0;
    //V_T.ArgI, mwcca6c729_26bf_4d26_b929_1aa139452cc5, index: 35
    //Unit: mole^(1)metre^(-3)
    _species_var[35] = PFILE(49) * 1000000.0;
    //V_T.ENT, mw407aec6d_b728_4aec_a53e_b75ab77d05f1, index: 36
    //Unit: mole^(1)metre^(-3)
    _species_var[36] = PFILE(50) * 1000000.0;
    //V_LN.nT0, mwebc60bab_74db_4cea_b897_335ba47a6603, index: 37
    //Unit: mole^(1)
    _species_var[37] = PFILE(51) * 1.66053872801495e-24;
    //V_LN.aT0, mw191547e0_918d_40cc_8fd8_9b683bb54293, index: 38
    //Unit: mole^(1)
    _species_var[38] = PFILE(52) * 1.66053872801495e-24;
    //V_LN.T0, mwf44705b4_092e_44fc_85fc_cf33ecad0d28, index: 39
    //Unit: mole^(1)
    _species_var[39] = PFILE(53) * 1.66053872801495e-24;
    //V_LN.IL2, mw6f1c494d_2927_4261_99e2_28f145f431c9, index: 40
    //Unit: mole^(1)metre^(-3)
    _species_var[40] = PFILE(54) * 999999999.9999999;
    //V_LN.nT1, mwcf120fb7_0510_49c1_96c4_f69eec93f4e1, index: 41
    //Unit: mole^(1)
    _species_var[41] = PFILE(55) * 1.66053872801495e-24;
    //V_LN.aT1, mw6ff1f892_a9b4_46bc_acdf_1a640e96eac7, index: 42
    //Unit: mole^(1)
    _species_var[42] = PFILE(56) * 1.66053872801495e-24;
    //V_LN.T1, mw58d7bb5f_2358_4c5b_bf9d_18dd7703a714, index: 43
    //Unit: mole^(1)
    _species_var[43] = PFILE(57) * 1.66053872801495e-24;
    //V_LN.APC, mw5a98782c_dfde_4296_8212_90500fa83148, index: 44
    //Unit: mole^(1)
    _species_var[44] = PFILE(58) * 1.66053872801495e-24;
    //V_LN.mAPC, mw6e2081f5_2ee1_4152_9153_d8b7f7373b50, index: 45
    //Unit: mole^(1)
    _species_var[45] = PFILE(59) * 1.66053872801495e-24;
    //V_LN.P0, mw73704559_1a6a_4e49_90c7_e8649c1d73b7, index: 46
    //Unit: mole^(1)metre^(-3)
    _species_var[46] = PFILE(60) * 999999999.9999999;
    //V_LN.P1, mwb7c54f02_3a02_497c_a2bd_6bb01949a193, index: 47
    //Unit: mole^(1)metre^(-3)
    _species_var[47] = PFILE(61) * 999999999.9999999;
    //V_LN.nivo, mw2369c64a_96ba_4362_86a1_d63fe83c8e42, index: 48
    //Unit: mole^(1)metre^(-3)
    _species_var[48] = PFILE(62) * 999999999.9999999;
    //V_LN.durv, mw6e060274_2daf_449c_a80b_44d9b16106c2, index: 49
    //Unit: mole^(1)metre^(-3)
    _species_var[49] = PFILE(63) * 999999999.9999999;
    //V_LN.ipi, mwd80c7730_4ed6_48b2_bb2f_3e748081d682, index: 50
    //Unit: mole^(1)metre^(-3)
    _species_var[50] = PFILE(64) * 999999999.9999999;
    //V_LN.ENT, mwa4006c5c_a807_4568_a3e6_d57e4a63d7fb, index: 51
    //Unit: mole^(1)metre^(-3)
    _species_var[51] = PFILE(65) * 999999999.9999999;
    //V_e.P0, mw6e212437_7abd_46d3_863a_5e5a22368719, index: 52
    //Unit: mole^(1)metre^(-3)
    _species_var[52] = PFILE(66) * 999.9999999999999;
    //V_e.p0, mwd341a5c2_449a_468b_a547_0198fcbe0903, index: 53
    //Unit: mole^(1)metre^(-3)
    _species_var[53] = PFILE(67) * 999.9999999999999;
    //V_e.P1, mwc0f4f373_9723_40e8_b022_973e5611a533, index: 54
    //Unit: mole^(1)metre^(-3)
    _species_var[54] = PFILE(68) * 999.9999999999999;
    //V_e.p1, mwd101a415_049d_4e9a_8c69_d54d6768a442, index: 55
    //Unit: mole^(1)metre^(-3)
    _species_var[55] = PFILE(69) * 999.9999999999999;
    //A_e.M1, mwc3ece17c_334f_4435_8944_477e8e943121, index: 56
    //Unit: mole^(1)metre^(-2)
    _species_var[56] = PFILE(70) * 1.66053872801495e-12;
    //A_e.M1p0, mw16bec048_b1a5_4c69_b34f_27d405f542f0, index: 57
    //Unit: mole^(1)metre^(-2)
    _species_var[57] = PFILE(71) * 1.66053872801495e-12;
    //A_e.M1p1, mwa0c42763_5110_4f37_bd22_ee490d8aa032, index: 58
    //Unit: mole^(1)metre^(-2)
    _species_var[58] = PFILE(72) * 1.66053872801495e-12;
    //A_s.M1, mwc2c4ce95_bb5f_45d8_946e_a2d4dee68c12, index: 59
    //Unit: mole^(1)metre^(-2)
    _species_var[59] = PFILE(73) * 1.66053872801495e-12;
    //A_s.M1p0, mw55e72a69_3002_4096_870f_76e134d684a9, index: 60
    //Unit: mole^(1)metre^(-2)
    _species_var[60] = PFILE(74) * 1.66053872801495e-12;
    //A_s.M1p1, mwc7771359_cc99_4940_9323_d990493fe40d, index: 61
    //Unit: mole^(1)metre^(-2)
    _species_var[61] = PFILE(75) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL1, mwa285b215_0bbb_4d93_bd49_45559fbf778d, index: 62
    //Unit: mole^(1)metre^(-2)
    _species_var[62] = PFILE(76) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL2, mw2ff9df15_19f4_4956_8c90_696dc1c42b87, index: 63
    //Unit: mole^(1)metre^(-2)
    _species_var[63] = PFILE(77) * 1.66053872801495e-12;
    //syn_T_C1.PD1, mwf52b6036_afdb_4e09_a9fe_7a2f3b7af98c, index: 64
    //Unit: mole^(1)metre^(-2)
    _species_var[64] = PFILE(78) * 1.66053872801495e-12;
    //syn_T_C1.PDL1, mw87e73b5e_4f4a_4e1b_aad7_30c4e59228c5, index: 65
    //Unit: mole^(1)metre^(-2)
    _species_var[65] = PFILE(79) * 1.66053872801495e-12;
    //syn_T_C1.PDL2, mw36656e4d_d5d2_4add_b70b_574adecafe9b, index: 66
    //Unit: mole^(1)metre^(-2)
    _species_var[66] = PFILE(80) * 1.66053872801495e-12;
    //syn_T_C1.PD1_nivo, mwcdf36a32_52dc_47ed_922e_22b59399f2bf, index: 67
    //Unit: mole^(1)metre^(-2)
    _species_var[67] = PFILE(81) * 1.66053872801495e-12;
    //syn_T_C1.PD1_nivo_PD1, mw2c2e61ba_da7b_47ef_a389_b004b272fcd8, index: 68
    //Unit: mole^(1)metre^(-2)
    _species_var[68] = PFILE(82) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_durv, mw8dba18a6_e88c_4c9d_ae32_8218f2427afc, index: 69
    //Unit: mole^(1)metre^(-2)
    _species_var[69] = PFILE(83) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_durv_PDL1, mwb12b733c_212a_4352_9c20_f3e6365f7356, index: 70
    //Unit: mole^(1)metre^(-2)
    _species_var[70] = PFILE(84) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1, mwf74be379_9b70_4aeb_88d5_8dc63cb481a0, index: 71
    //Unit: mole^(1)metre^(-2)
    _species_var[71] = PFILE(85) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_durv, mw927433f0_22fe_4de1_b2e7_3988f53dafd4, index: 72
    //Unit: mole^(1)metre^(-2)
    _species_var[72] = PFILE(86) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_durv_TPDL1, mw08b60496_630a_4840_bf49_028afa3203db, index: 73
    //Unit: mole^(1)metre^(-2)
    _species_var[73] = PFILE(87) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80, mw47782763_66ef_4a9b_bfe3_0e7e0873913e, index: 74
    //Unit: mole^(1)metre^(-2)
    _species_var[74] = PFILE(88) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80_CD28, mw74f8e7f0_b3a9_42af_9164_3fc4686247a0, index: 75
    //Unit: mole^(1)metre^(-2)
    _species_var[75] = PFILE(89) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD86, mw46989c2c_5a12_474f_987a_3e480a442f78, index: 76
    //Unit: mole^(1)metre^(-2)
    _species_var[76] = PFILE(90) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4, mw8885eb49_c385_444f_8b0a_6fe53ee7283f, index: 77
    //Unit: mole^(1)metre^(-2)
    _species_var[77] = PFILE(91) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80, mw6322ed31_01e3_4e93_af12_f35af304a5d2, index: 78
    //Unit: mole^(1)metre^(-2)
    _species_var[78] = PFILE(92) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_CD80_CTLA4, mw3254bf69_8762_4c12_935e_6bb9dabd1491, index: 79
    //Unit: mole^(1)metre^(-2)
    _species_var[79] = PFILE(93) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80_CTLA4, mw5952b0f0_ea30_47e2_87a9_0522377b13bd, index: 80
    //Unit: mole^(1)metre^(-2)
    _species_var[80] = PFILE(94) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4, mw0a4c4a3d_d2b1_4a3e_a0d7_5f6fb1eb111d, index: 81
    //Unit: mole^(1)metre^(-2)
    _species_var[81] = PFILE(95) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4_CD86, mw76d61402_92cc_4953_93d9_124a45b5dc93, index: 82
    //Unit: mole^(1)metre^(-2)
    _species_var[82] = PFILE(96) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_CD80, mwdd30eae0_06ac_43b0_878f_3897ebfee6c5, index: 83
    //Unit: mole^(1)metre^(-2)
    _species_var[83] = PFILE(97) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_CD80_TPDL1, mwbc4bab2f_4f13_40c5_92fb_e8a72445599f, index: 84
    //Unit: mole^(1)metre^(-2)
    _species_var[84] = PFILE(98) * 1.66053872801495e-12;
    //syn_T_C1.CD28, mw704d6065_fad3_49a5_a0e6_a928bff19ee1, index: 85
    //Unit: mole^(1)metre^(-2)
    _species_var[85] = PFILE(99) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4, mw1f32efee_79a0_486a_8a7a_1bf653a1546b, index: 86
    //Unit: mole^(1)metre^(-2)
    _species_var[86] = PFILE(100) * 1.66053872801495e-12;
    //syn_T_C1.CD80, mw43e8de59_dcb9_45ce_9281_dec0db54165e, index: 87
    //Unit: mole^(1)metre^(-2)
    _species_var[87] = PFILE(101) * 1.66053872801495e-12;
    //syn_T_C1.CD86, mw9d2ab5be_0107_4617_929f_fd9e0e57fa0d, index: 88
    //Unit: mole^(1)metre^(-2)
    _species_var[88] = PFILE(102) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_ipi, mw69f6a218_1e87_4157_ace0_116b2ae14d96, index: 89
    //Unit: mole^(1)metre^(-2)
    _species_var[89] = PFILE(103) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_ipi_CTLA4, mwcabd5c26_6313_4411_a44e_76cae7a474d9, index: 90
    //Unit: mole^(1)metre^(-2)
    _species_var[90] = PFILE(104) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL1, mw47ef397c_30e3_4eec_9230_ee36f3278bb4, index: 91
    //Unit: mole^(1)metre^(-2)
    _species_var[91] = PFILE(105) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL2, mwfb7f4630_b154_4fae_a6d6_d71461fdf4f2, index: 92
    //Unit: mole^(1)metre^(-2)
    _species_var[92] = PFILE(106) * 1.66053872801495e-12;
    //syn_T_APC.PD1, mw13a74242_41f3_4e62_b7b1_a2ccc3a3ac94, index: 93
    //Unit: mole^(1)metre^(-2)
    _species_var[93] = PFILE(107) * 1.66053872801495e-12;
    //syn_T_APC.PDL1, mw37d3fe0c_0c28_4071_ab08_edd43d7de4c4, index: 94
    //Unit: mole^(1)metre^(-2)
    _species_var[94] = PFILE(108) * 1.66053872801495e-12;
    //syn_T_APC.PDL2, mwdaa6097a_1baf_4f1a_a0b0_aae346c8b612, index: 95
    //Unit: mole^(1)metre^(-2)
    _species_var[95] = PFILE(109) * 1.66053872801495e-12;
    //syn_T_APC.PD1_nivo, mw122d8fbc_6fa1_40b3_b9c2_35b7dc70cd3f, index: 96
    //Unit: mole^(1)metre^(-2)
    _species_var[96] = PFILE(110) * 1.66053872801495e-12;
    //syn_T_APC.PD1_nivo_PD1, mw28be4822_13e3_4c26_8e75_359bda2b49a8, index: 97
    //Unit: mole^(1)metre^(-2)
    _species_var[97] = PFILE(111) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_durv, mw71e39d87_b022_4538_9d84_b6880d6cd4c0, index: 98
    //Unit: mole^(1)metre^(-2)
    _species_var[98] = PFILE(112) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_durv_PDL1, mw6dd02de7_3ac2_473b_a571_99d8f504a2a6, index: 99
    //Unit: mole^(1)metre^(-2)
    _species_var[99] = PFILE(113) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1, mw818d0b31_c40e_4340_944a_75c058e64e08, index: 100
    //Unit: mole^(1)metre^(-2)
    _species_var[100] = PFILE(114) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_durv, mw40481391_a71d_4ffe_b6a3_95b2104ea094, index: 101
    //Unit: mole^(1)metre^(-2)
    _species_var[101] = PFILE(115) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_durv_TPDL1, mw8369aaf6_5bcc_4866_97b6_6ef8e4d3c91b, index: 102
    //Unit: mole^(1)metre^(-2)
    _species_var[102] = PFILE(116) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80, mwc79e7416_6ad5_4bb1_8476_ed5d17ad53ce, index: 103
    //Unit: mole^(1)metre^(-2)
    _species_var[103] = PFILE(117) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80_CD28, mw4547854f_e867_430b_b343_97453aa59678, index: 104
    //Unit: mole^(1)metre^(-2)
    _species_var[104] = PFILE(118) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD86, mwcbb44837_205a_40f2_bfd6_d9981dcb2c94, index: 105
    //Unit: mole^(1)metre^(-2)
    _species_var[105] = PFILE(119) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4, mw40cfc0ca_4762_4a89_838a_66deb1eacd70, index: 106
    //Unit: mole^(1)metre^(-2)
    _species_var[106] = PFILE(120) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80, mw8ad55d8a_d6df_4a43_adca_2e74efb19f37, index: 107
    //Unit: mole^(1)metre^(-2)
    _species_var[107] = PFILE(121) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_CD80_CTLA4, mw19e00d1d_bf7e_4d3d_a104_3f496d8c6cdf, index: 108
    //Unit: mole^(1)metre^(-2)
    _species_var[108] = PFILE(122) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80_CTLA4, mw5473d67e_62c2_4d15_9aec_2192a55de156, index: 109
    //Unit: mole^(1)metre^(-2)
    _species_var[109] = PFILE(123) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4, mwba275264_896a_42ed_bfd0_b9d438080161, index: 110
    //Unit: mole^(1)metre^(-2)
    _species_var[110] = PFILE(124) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4_CD86, mwa8d37cce_ee58_4e86_a1ba_09bd44129532, index: 111
    //Unit: mole^(1)metre^(-2)
    _species_var[111] = PFILE(125) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_CD80, mw263da3d6_5007_4039_b8da_05ec1d6d1de4, index: 112
    //Unit: mole^(1)metre^(-2)
    _species_var[112] = PFILE(126) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_CD80_TPDL1, mw8e409b8f_51a1_4717_aa8e_09d1cea652f5, index: 113
    //Unit: mole^(1)metre^(-2)
    _species_var[113] = PFILE(127) * 1.66053872801495e-12;
    //syn_T_APC.CD28, mw8dddbd87_7cba_4480_aacf_322fefe15e00, index: 114
    //Unit: mole^(1)metre^(-2)
    _species_var[114] = PFILE(128) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4, mwad31dd04_8f1c_468b_bf95_1024c36e69d2, index: 115
    //Unit: mole^(1)metre^(-2)
    _species_var[115] = PFILE(129) * 1.66053872801495e-12;
    //syn_T_APC.CD80, mwbfb8df59_75b6_4449_b6a2_bb5133f35c35, index: 116
    //Unit: mole^(1)metre^(-2)
    _species_var[116] = PFILE(130) * 1.66053872801495e-12;
    //syn_T_APC.CD86, mw2f62e6c4_a8f9_4511_b779_ef2dfc381412, index: 117
    //Unit: mole^(1)metre^(-2)
    _species_var[117] = PFILE(131) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_ipi, mw48ad0e90_4bc8_4c2d_b3b1_bce9f249e72d, index: 118
    //Unit: mole^(1)metre^(-2)
    _species_var[118] = PFILE(132) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_ipi_CTLA4, mwb86eb755_5b57_4338_b9db_0d5da84f0b05, index: 119
    //Unit: mole^(1)metre^(-2)
    _species_var[119] = PFILE(133) * 1.66053872801495e-12;
    
    return;
}

void ODE_system::setup_instance_tolerance(Param& param){

    //Tolerance
    realtype reltol = PFILE(3);
    realtype abstol_base = PFILE(4);
    N_Vector abstol = N_VNew_Serial(_neq);

    for (size_t i = 0; i < 120; i++)
    {
        NV_DATA_S(abstol)[i] = abstol_base * get_unit_conversion_species(i);
    }
    int flag = CVodeSVtolerances(_cvode_mem, reltol, abstol);
    check_flag(&flag, "CVodeSVtolerances", 1);

    
    return;
}

void ODE_system::eval_init_assignment(void){
    //Assignment Rules required before IA
    //InitialAssignment
    _species_var[64] = _class_parameter[146] / _class_parameter[72];
    _species_var[85] = _class_parameter[147] / _class_parameter[72];
    _species_var[86] = _class_parameter[148] / _class_parameter[71];
    _species_var[71] = _class_parameter[149] / _class_parameter[72];
    _species_var[65] = _class_parameter[150] / _class_parameter[73];
    _species_var[66] = _class_parameter[151] / _class_parameter[73];
    _species_var[87] = _class_parameter[152] / _class_parameter[73];
    _species_var[88] = _class_parameter[153] / _class_parameter[73];
    _species_var[93] = _class_parameter[146] / _class_parameter[72];
    _species_var[114] = _class_parameter[147] / _class_parameter[72];
    _species_var[115] = _class_parameter[148] / _class_parameter[71];
    _species_var[100] = _class_parameter[149] / _class_parameter[72];
    _species_var[94] = _class_parameter[154] / _class_parameter[74];
    _species_var[95] = _class_parameter[155] / _class_parameter[74];
    _species_var[116] = _class_parameter[156] / _class_parameter[74];
    _species_var[117] = _class_parameter[157] / _class_parameter[74];
    _species_var[29] = _class_parameter[158];
    _species_var[14] = _class_parameter[158];

    updateVar();
    
    return;
}
void ODE_system::setupEvents(void){

    _nevent = 1;
    _nroot = 1;

    _trigger_element_type = std::vector<EVENT_TRIGGER_ELEM_TYPE>(_nroot, TRIGGER_NON_INSTANT);
    _trigger_element_satisfied = std::vector<bool>(_nroot, false);
    _event_triggered = std::vector<bool>(_nevent, false);

    //V_T.C1 < (0.9 * cell)
    _trigger_element_type[0] = TRIGGER_NON_INSTANT;

    _event_triggered[0] = true;

    return;
}
int ODE_system::f(realtype t, N_Vector y, N_Vector ydot, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    realtype AUX_VAR_V_T = PARAM(13) + PARAM(11) * SPVAR(18) + PARAM(12) * SPVAR(19) + PARAM(11) * SPVAR(20) + PARAM(12) * SPVAR(21) + PARAM(12) * SPVAR(22);

    realtype AUX_VAR_C_total = 0.0 * PARAM(9) + SPVAR(20);

    realtype AUX_VAR_T_total = 0.0 * PARAM(9) + SPVAR(21) + SPVAR(22);

    realtype AUX_VAR_T_total_LN = 0.0 * PARAM(9) + SPVAR(43);

    realtype AUX_VAR_H_PD1_C1 = std::pow((SPVAR(62) + SPVAR(63)) / PARAM(142), PARAM(143)) / (std::pow((SPVAR(62) + SPVAR(63)) / PARAM(142), PARAM(143)) + 1.0);

    realtype AUX_VAR_H_MDSC_C1 = 1.0 - (1.0 - SPVAR(34) / (PARAM(174) + SPVAR(34))) * (1.0 - SPVAR(35) / (PARAM(173) + SPVAR(35))) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(181)));

    realtype AUX_VAR_Tregs_ = SPVAR(21);

    realtype AUX_VAR_H_CD28_APC = std::pow((SPVAR(103) + SPVAR(105) + 2.0 * SPVAR(104)) / PARAM(144), PARAM(145)) / (std::pow((SPVAR(103) + SPVAR(105) + 2.0 * SPVAR(104)) / PARAM(144), PARAM(145)) + 1.0);

    realtype AUX_VAR_pTCR_p0_MHC_tot = 0.5 * (SPVAR(60) / PARAM(19) + PARAM(77) + PARAM(76) / PARAM(75) - PARAM(77) * std::pow(std::pow((SPVAR(60) / PARAM(19) + PARAM(77) + PARAM(76) / PARAM(75)) / PARAM(77), 2.0) - 4.0 * SPVAR(60) / PARAM(19) / PARAM(77), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p1_MHC_tot = PARAM(87) / (PARAM(87) + PARAM(89)) * std::pow(PARAM(88) / (PARAM(87) + PARAM(88)), PARAM(90)) * 0.5 * (SPVAR(61) / PARAM(38) + PARAM(91) + PARAM(87) / PARAM(86) - PARAM(91) * std::pow(std::pow((SPVAR(61) / PARAM(38) + PARAM(91) + PARAM(87) / PARAM(86)) / PARAM(91), 2.0) - 4.0 * SPVAR(61) / PARAM(38) / PARAM(91), 1.0 / 2.0));

    realtype AUX_VAR_H_CD28_C1 = std::pow((SPVAR(74) + SPVAR(76) + 2.0 * SPVAR(75)) / PARAM(144), PARAM(145)) / (std::pow((SPVAR(74) + SPVAR(76) + 2.0 * SPVAR(75)) / PARAM(144), PARAM(145)) + 1.0);

    realtype AUX_VAR_H_PD1_APC = std::pow((SPVAR(91) + SPVAR(92)) / PARAM(142), PARAM(143)) / (std::pow((SPVAR(91) + SPVAR(92)) / PARAM(142), PARAM(143)) + 1.0);

    realtype AUX_VAR_H_Treg_T = std::pow((SPVAR(30) + 2.0 * SPVAR(31)) / PARAM(159), PARAM(160)) / (std::pow((SPVAR(30) + 2.0 * SPVAR(31)) / PARAM(159), PARAM(160)) + 1.0);

    realtype AUX_VAR_H_Treg_P = std::pow((SPVAR(15) + 2.0 * SPVAR(16)) / PARAM(159), PARAM(160)) / (std::pow((SPVAR(15) + 2.0 * SPVAR(16)) / PARAM(159), PARAM(160)) + 1.0);

    realtype AUX_VAR_H_ENT_C1 = SPVAR(36) / (SPVAR(36) + PARAM(164));

    realtype AUX_VAR_H_APC = PARAM(60) * SPVAR(44) / (PARAM(60) * SPVAR(44) + AUX_VAR_T_total_LN + PARAM(9));

    realtype AUX_VAR_H_mAPC = PARAM(60) * SPVAR(45) / (PARAM(60) * SPVAR(45) + AUX_VAR_T_total_LN + PARAM(9));

    realtype AUX_VAR_R_Tcell = 0.0 * PARAM(9) / PARAM(10) + PARAM(49) * SPVAR(22) * SPVAR(20) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype AUX_VAR_N_aT = PARAM(34) + PARAM(35) * AUX_VAR_H_CD28_APC + PARAM(36) * SPVAR(40) / (PARAM(32) + SPVAR(40));

    realtype AUX_VAR_H_P0 = AUX_VAR_pTCR_p0_MHC_tot / (AUX_VAR_pTCR_p0_MHC_tot + PARAM(69));

    realtype AUX_VAR_H_P1 = AUX_VAR_pTCR_p1_MHC_tot / (AUX_VAR_pTCR_p1_MHC_tot + PARAM(84));

    //Reaction fluxes:

    realtype ReactionFlux1 = PARAM(8) * SPVAR(18);

    realtype ReactionFlux2 = PARAM(8) * SPVAR(19);

    realtype ReactionFlux3 = PARAM(14) * SPVAR(20) * (1.0 - AUX_VAR_C_total / PARAM(15)) * (1.0 - AUX_VAR_H_ENT_C1);

    realtype ReactionFlux4 = PARAM(16) * SPVAR(20);

    realtype ReactionFlux5 = PARAM(20) * PARAM(19);

    realtype ReactionFlux6 = PARAM(21) * SPVAR(37);

    realtype ReactionFlux7 = PARAM(22) * AUX_VAR_H_APC * AUX_VAR_H_P0 * SPVAR(37);

    realtype ReactionFlux8 = PARAM(23) / AUX_VAR_N_aT * SPVAR(38);

    realtype ReactionFlux9 = PARAM(23) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(38);

    realtype ReactionFlux10 = PARAM(24) * SPVAR(0);

    realtype ReactionFlux11 = PARAM(24) * SPVAR(9);

    realtype ReactionFlux12 = PARAM(24) * SPVAR(21);

    realtype ReactionFlux13 = PARAM(24) * SPVAR(39);

    realtype ReactionFlux14 = PARAM(25) * SPVAR(0);

    realtype ReactionFlux15 = PARAM(26) * SPVAR(9);

    realtype ReactionFlux16 = PARAM(27) * AUX_VAR_V_T * SPVAR(0) * (1.0 - SPVAR(21) / (PARAM(178) * AUX_VAR_V_T));      

    realtype ReactionFlux17 = PARAM(28) * SPVAR(39);

    realtype ReactionFlux18 = PARAM(29) * SPVAR(40) * PARAM(2);

    realtype ReactionFlux19 = PARAM(30) * AUX_VAR_T_total_LN * SPVAR(40) / (PARAM(32) + SPVAR(40));

    realtype ReactionFlux20 = PARAM(30) * SPVAR(39) * SPVAR(40) / (PARAM(33) + SPVAR(40));

    realtype ReactionFlux21 = PARAM(39) * PARAM(38);

    realtype ReactionFlux22 = PARAM(40) * SPVAR(41);

    realtype ReactionFlux23 = PARAM(41) * AUX_VAR_H_mAPC * AUX_VAR_H_P1 * SPVAR(41);

    realtype ReactionFlux24 = PARAM(42) / AUX_VAR_N_aT * SPVAR(42);

    realtype ReactionFlux25 = PARAM(42) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(42);

    realtype ReactionFlux26 = PARAM(43) * SPVAR(1);

    realtype ReactionFlux27 = PARAM(43) * SPVAR(10);

    realtype ReactionFlux28 = PARAM(43) * SPVAR(22);

    realtype ReactionFlux29 = PARAM(43) * SPVAR(43);

    realtype ReactionFlux30 = PARAM(37) * SPVAR(22) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux31 = PARAM(48) * SPVAR(22) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux32 = PARAM(44) * SPVAR(1);

    realtype ReactionFlux33 = PARAM(45) * SPVAR(10);

    realtype ReactionFlux34 = PARAM(46) * AUX_VAR_V_T * SPVAR(1);

    realtype ReactionFlux35 = PARAM(47) * SPVAR(43);

    realtype ReactionFlux36 = PARAM(31) * SPVAR(42);

    realtype ReactionFlux37 = PARAM(49) * SPVAR(22) * SPVAR(20) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux38 = PARAM(52) * (PARAM(54) * AUX_VAR_V_T - SPVAR(23));

    realtype ReactionFlux39 = PARAM(52) * (PARAM(55) * PARAM(2) - SPVAR(44));

    realtype ReactionFlux40 = PARAM(50) * SPVAR(25) / (SPVAR(25) + PARAM(58)) * SPVAR(23);

    realtype ReactionFlux41 = PARAM(51) * SPVAR(24);

    realtype ReactionFlux42 = PARAM(53) * SPVAR(24);

    realtype ReactionFlux43 = PARAM(53) * SPVAR(45);

    realtype ReactionFlux44 = PARAM(56) * (PARAM(57) - SPVAR(25)) * AUX_VAR_V_T;

    realtype ReactionFlux45 = AUX_VAR_R_Tcell * PARAM(59);

    realtype ReactionFlux46 = PARAM(62) * SPVAR(56) * PARAM(4) - PARAM(61) * SPVAR(59) * PARAM(5);

    realtype ReactionFlux47 = PARAM(19) * PARAM(70) * (PARAM(16) + PARAM(17) + PARAM(49) * SPVAR(22) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(20) * AUX_VAR_V_T;

    realtype ReactionFlux48 = PARAM(64) * SPVAR(46) * PARAM(2);

    realtype ReactionFlux49 = PARAM(63) * SPVAR(45) * SPVAR(46) * PARAM(2);

    realtype ReactionFlux50 = PARAM(63) * PARAM(9) * SPVAR(46) * PARAM(3);

    realtype ReactionFlux51 = PARAM(65) * SPVAR(52) * PARAM(3);

    realtype ReactionFlux52 = PARAM(66) * SPVAR(53) * PARAM(3);

    realtype ReactionFlux53 = PARAM(67) * SPVAR(53) * SPVAR(56) * PARAM(4);

    realtype ReactionFlux54 = PARAM(68) * PARAM(67) * SPVAR(57) * PARAM(4);

    realtype ReactionFlux55 = PARAM(68) * PARAM(67) * SPVAR(60) * PARAM(5);

    realtype ReactionFlux56 = PARAM(62) * SPVAR(57) * PARAM(4);

    realtype ReactionFlux57 = PARAM(38) * PARAM(85) * (PARAM(16) + PARAM(17) + PARAM(49) * SPVAR(22) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(20) * AUX_VAR_V_T;

    realtype ReactionFlux58 = PARAM(79) * SPVAR(47) * PARAM(2);

    realtype ReactionFlux59 = PARAM(78) * SPVAR(45) * SPVAR(47) * PARAM(2);

    realtype ReactionFlux60 = PARAM(78) * PARAM(9) * SPVAR(47) * PARAM(3);

    realtype ReactionFlux61 = PARAM(80) * SPVAR(54) * PARAM(3);

    realtype ReactionFlux62 = PARAM(81) * SPVAR(55) * PARAM(3);

    realtype ReactionFlux63 = PARAM(82) * SPVAR(55) * SPVAR(56) * PARAM(4);

    realtype ReactionFlux64 = PARAM(83) * PARAM(82) * SPVAR(58) * PARAM(4);

    realtype ReactionFlux65 = PARAM(83) * PARAM(82) * SPVAR(61) * PARAM(5);

    realtype ReactionFlux66 = PARAM(62) * SPVAR(58) * PARAM(4);

    realtype ReactionFlux67 = PARAM(93) * (SPVAR(2) / PARAM(98) - SPVAR(11) / PARAM(99));

    realtype ReactionFlux68 = PARAM(94) * (SPVAR(2) / PARAM(98) - SPVAR(26) / PARAM(100));

    realtype ReactionFlux69 = PARAM(95) * (SPVAR(2) / PARAM(98) - SPVAR(48) / PARAM(101));

    realtype ReactionFlux70 = PARAM(96) * SPVAR(26) / PARAM(100) * AUX_VAR_V_T;

    realtype ReactionFlux71 = PARAM(96) * SPVAR(48) / PARAM(101) * PARAM(2);

    realtype ReactionFlux72 = PARAM(97) * SPVAR(2);

    realtype ReactionFlux73 = PARAM(102) * (SPVAR(3) / PARAM(107) - SPVAR(12) / PARAM(108));

    realtype ReactionFlux74 = PARAM(103) * (SPVAR(3) / PARAM(107) - SPVAR(27) / PARAM(109));

    realtype ReactionFlux75 = PARAM(104) * (SPVAR(3) / PARAM(107) - SPVAR(49) / PARAM(110));

    realtype ReactionFlux76 = PARAM(105) * SPVAR(27) / PARAM(109) * AUX_VAR_V_T;

    realtype ReactionFlux77 = PARAM(105) * SPVAR(49) / PARAM(110) * PARAM(2);

    realtype ReactionFlux78 = PARAM(106) * SPVAR(3);

    realtype ReactionFlux79 = PARAM(111) * (SPVAR(4) / PARAM(116) - SPVAR(13) / PARAM(117));

    realtype ReactionFlux80 = PARAM(112) * (SPVAR(4) / PARAM(116) - SPVAR(28) / PARAM(118));

    realtype ReactionFlux81 = PARAM(113) * (SPVAR(4) / PARAM(116) - SPVAR(50) / PARAM(119));

    realtype ReactionFlux82 = PARAM(114) * SPVAR(28) / PARAM(118) * AUX_VAR_V_T;

    realtype ReactionFlux83 = PARAM(114) * SPVAR(50) / PARAM(119) * PARAM(2);

    realtype ReactionFlux84 = PARAM(115) * SPVAR(4);

    realtype ReactionFlux85 = (PARAM(92) * SPVAR(64) * SPVAR(65) - PARAM(129) * SPVAR(62)) * PARAM(6);

    realtype ReactionFlux86 = (PARAM(120) * SPVAR(64) * SPVAR(66) - PARAM(130) * SPVAR(63)) * PARAM(6);

    realtype ReactionFlux87 = (2.0 * PARAM(121) * (SPVAR(64) * SPVAR(26) / PARAM(100)) - PARAM(131) * SPVAR(67)) * PARAM(6);

    realtype ReactionFlux88 = (PARAM(139) * PARAM(121) * SPVAR(64) * SPVAR(67) - 2.0 * PARAM(131) * SPVAR(68)) * PARAM(6);

    realtype ReactionFlux89 = (2.0 * PARAM(122) * (SPVAR(65) * SPVAR(27) / PARAM(109)) - PARAM(132) * SPVAR(69)) * PARAM(6);

    realtype ReactionFlux90 = (PARAM(140) * PARAM(122) * SPVAR(65) * SPVAR(69) - 2.0 * PARAM(132) * SPVAR(70)) * PARAM(6);

    realtype ReactionFlux91 = (2.0 * PARAM(123) * SPVAR(85) * SPVAR(87) - PARAM(133) * SPVAR(74)) * PARAM(6);

    realtype ReactionFlux92 = (PARAM(123) * SPVAR(85) * SPVAR(74) - 2.0 * PARAM(133) * SPVAR(75)) * PARAM(6);

    realtype ReactionFlux93 = (PARAM(124) * SPVAR(85) * SPVAR(88) - PARAM(134) * SPVAR(76)) * PARAM(6);

    realtype ReactionFlux94 = (4.0 * PARAM(125) * SPVAR(86) * SPVAR(87) - PARAM(135) * SPVAR(77)) * PARAM(6);

    realtype ReactionFlux95 = (PARAM(125) * SPVAR(86) * SPVAR(77) - 2.0 * PARAM(135) * SPVAR(79)) * PARAM(6);

    realtype ReactionFlux96 = (PARAM(125) * SPVAR(87) * SPVAR(79) - PARAM(135) * SPVAR(80)) * PARAM(6);

    realtype ReactionFlux97 = (PARAM(125) * SPVAR(77) * SPVAR(87) - 2.0 * PARAM(135) * SPVAR(78)) * PARAM(6);

    realtype ReactionFlux98 = (PARAM(125) * SPVAR(86) * SPVAR(78) - PARAM(135) * SPVAR(80)) * PARAM(6);

    realtype ReactionFlux99 = (2.0 * PARAM(126) * SPVAR(86) * SPVAR(88) - PARAM(136) * SPVAR(81)) * PARAM(6);

    realtype ReactionFlux100 = (PARAM(126) * SPVAR(81) * SPVAR(88) - 2.0 * PARAM(136) * SPVAR(82)) * PARAM(6);

    realtype ReactionFlux101 = (4.0 * PARAM(128) * (SPVAR(86) * SPVAR(28) / PARAM(118)) - PARAM(138) * SPVAR(89)) * PARAM(6);

    realtype ReactionFlux102 = (PARAM(141) * PARAM(128) * SPVAR(86) * SPVAR(89) - 2.0 * PARAM(138) * SPVAR(90)) * PARAM(6);

    realtype ReactionFlux103 = (2.0 * PARAM(127) * SPVAR(87) * SPVAR(71) - PARAM(137) * SPVAR(83)) * PARAM(6);

    realtype ReactionFlux104 = (PARAM(127) * SPVAR(83) * SPVAR(71) - 2.0 * PARAM(137) * SPVAR(84)) * PARAM(6);

    realtype ReactionFlux105 = (2.0 * PARAM(122) * (SPVAR(71) * SPVAR(27) / PARAM(109)) - PARAM(132) * SPVAR(72)) * PARAM(6);

    realtype ReactionFlux106 = (PARAM(140) * PARAM(122) * SPVAR(71) * SPVAR(72) - 2.0 * PARAM(132) * SPVAR(73)) * PARAM(6);

    realtype ReactionFlux107 = (PARAM(92) * SPVAR(93) * SPVAR(94) - PARAM(129) * SPVAR(91)) * PARAM(7);

    realtype ReactionFlux108 = (PARAM(120) * SPVAR(93) * SPVAR(95) - PARAM(130) * SPVAR(92)) * PARAM(7);

    realtype ReactionFlux109 = (2.0 * PARAM(121) * (SPVAR(93) * SPVAR(48) / PARAM(101)) - PARAM(131) * SPVAR(96)) * PARAM(7);

    realtype ReactionFlux110 = (PARAM(139) * PARAM(121) * SPVAR(93) * SPVAR(96) - 2.0 * PARAM(131) * SPVAR(97)) * PARAM(7);

    realtype ReactionFlux111 = (2.0 * PARAM(122) * (SPVAR(94) * SPVAR(49) / PARAM(110)) - PARAM(132) * SPVAR(98)) * PARAM(7);

    realtype ReactionFlux112 = (PARAM(140) * PARAM(122) * SPVAR(94) * SPVAR(98) - 2.0 * PARAM(132) * SPVAR(99)) * PARAM(7);

    realtype ReactionFlux113 = (2.0 * PARAM(123) * SPVAR(114) * SPVAR(116) - PARAM(133) * SPVAR(103)) * PARAM(7);

    realtype ReactionFlux114 = (PARAM(123) * SPVAR(114) * SPVAR(103) - 2.0 * PARAM(133) * SPVAR(104)) * PARAM(7);

    realtype ReactionFlux115 = (PARAM(124) * SPVAR(114) * SPVAR(117) - PARAM(134) * SPVAR(105)) * PARAM(7);

    realtype ReactionFlux116 = (4.0 * PARAM(125) * SPVAR(115) * SPVAR(116) - PARAM(135) * SPVAR(106)) * PARAM(7);

    realtype ReactionFlux117 = (PARAM(125) * SPVAR(115) * SPVAR(106) - 2.0 * PARAM(135) * SPVAR(108)) * PARAM(7);

    realtype ReactionFlux118 = (PARAM(125) * SPVAR(116) * SPVAR(108) - PARAM(135) * SPVAR(109)) * PARAM(7);

    realtype ReactionFlux119 = (PARAM(125) * SPVAR(106) * SPVAR(116) - 2.0 * PARAM(135) * SPVAR(107)) * PARAM(7);

    realtype ReactionFlux120 = (PARAM(125) * SPVAR(115) * SPVAR(107) - PARAM(135) * SPVAR(109)) * PARAM(7);

    realtype ReactionFlux121 = (2.0 * PARAM(126) * SPVAR(115) * SPVAR(117) - PARAM(136) * SPVAR(110)) * PARAM(7);

    realtype ReactionFlux122 = (PARAM(126) * SPVAR(110) * SPVAR(117) - 2.0 * PARAM(136) * SPVAR(111)) * PARAM(7);

    realtype ReactionFlux123 = (4.0 * PARAM(128) * (SPVAR(115) * SPVAR(50) / PARAM(119)) - PARAM(138) * SPVAR(118)) * PARAM(7);

    realtype ReactionFlux124 = (PARAM(141) * PARAM(128) * SPVAR(115) * SPVAR(118) - 2.0 * PARAM(138) * SPVAR(119)) * PARAM(7);

    realtype ReactionFlux125 = (2.0 * PARAM(127) * SPVAR(116) * SPVAR(100) - PARAM(137) * SPVAR(112)) * PARAM(7);

    realtype ReactionFlux126 = (PARAM(127) * SPVAR(112) * SPVAR(100) - 2.0 * PARAM(137) * SPVAR(113)) * PARAM(7);

    realtype ReactionFlux127 = (2.0 * PARAM(122) * (SPVAR(100) * SPVAR(49) / PARAM(110)) - PARAM(132) * SPVAR(101)) * PARAM(7);

    realtype ReactionFlux128 = (PARAM(140) * PARAM(122) * SPVAR(100) * SPVAR(101) - 2.0 * PARAM(132) * SPVAR(102)) * PARAM(7);

    realtype ReactionFlux129 = PARAM(128) * (SPVAR(29) * SPVAR(28) / PARAM(116)) - PARAM(138) * SPVAR(30);

    realtype ReactionFlux130 = PARAM(141) * PARAM(128) * SPVAR(29) * SPVAR(30) / PARAM(72) - PARAM(138) * SPVAR(31);

    realtype ReactionFlux131 = PARAM(128) * (SPVAR(14) * SPVAR(13) / PARAM(117)) - PARAM(138) * SPVAR(15);

    realtype ReactionFlux132 = PARAM(141) * PARAM(128) * SPVAR(14) * SPVAR(15) / PARAM(72) - PARAM(138) * SPVAR(16);

    realtype ReactionFlux133 = PARAM(161) * SPVAR(21) * AUX_VAR_H_Treg_T;

    realtype ReactionFlux134 = PARAM(161) * SPVAR(9) * AUX_VAR_H_Treg_P;

    realtype ReactionFlux135 = PARAM(162) * (PARAM(177) * AUX_VAR_V_T - SPVAR(32)) * (SPVAR(33) / (PARAM(175) + SPVAR(33)));

    realtype ReactionFlux136 = PARAM(180) * (PARAM(177) * AUX_VAR_V_T - SPVAR(32));      

    realtype ReactionFlux137 = PARAM(163) * SPVAR(32);

    realtype ReactionFlux138 = PARAM(165) * SPVAR(33) * AUX_VAR_V_T;

    realtype ReactionFlux139 = PARAM(166) * SPVAR(34) * AUX_VAR_V_T;

    realtype ReactionFlux140 = PARAM(167) * SPVAR(35) * AUX_VAR_V_T;

    realtype ReactionFlux141 = PARAM(168) * SPVAR(20) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(179)));

    realtype ReactionFlux142 = PARAM(169) * SPVAR(32) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(171)));

    realtype ReactionFlux143 = PARAM(170) * SPVAR(32) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(181)));

    realtype ReactionFlux144 = PARAM(172) * SPVAR(21) * (1.0 - SPVAR(21) / (PARAM(178) * AUX_VAR_V_T)) * SPVAR(35) / (PARAM(176) + SPVAR(35)) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(181)));

    realtype ReactionFlux145 = PARAM(188) * SPVAR(8) * PARAM(0);

    realtype ReactionFlux146 = PARAM(182) * SPVAR(6) * PARAM(0);

    realtype ReactionFlux147 = PARAM(183) * SPVAR(7) * PARAM(0);

    realtype ReactionFlux148 = PARAM(184) * SPVAR(5) / (SPVAR(5) + PARAM(185)) * PARAM(0);

    realtype ReactionFlux149 = PARAM(189) * (SPVAR(5) / PARAM(194) - SPVAR(17) / PARAM(195));

    realtype ReactionFlux150 = PARAM(190) * (SPVAR(5) / PARAM(194) - SPVAR(36) / PARAM(196));

    realtype ReactionFlux151 = PARAM(191) * (SPVAR(5) / PARAM(194) - SPVAR(51) / PARAM(197));

    realtype ReactionFlux152 = PARAM(192) * SPVAR(36) / PARAM(196) * AUX_VAR_V_T;

    realtype ReactionFlux153 = PARAM(192) * SPVAR(51) / PARAM(197) * PARAM(2);

    realtype ReactionFlux154 = PARAM(193) * SPVAR(5) * PARAM(0);
    
    //dydt:

    //d(V_C.T0)/dt

    if (SPVAR(21) / (PARAM(178) * AUX_VAR_V_T) < 1){
        realtype ReactionFlux16 = PARAM(27) * AUX_VAR_V_T * SPVAR(0) * (1.0 - SPVAR(21) / (PARAM(178) * AUX_VAR_V_T));      
        NV_DATA_S(ydot)[0] =  - ReactionFlux10 - ReactionFlux14 + ReactionFlux15 - ReactionFlux16 + ReactionFlux17;
    }
    else{
        realtype ReactionFlux16 = 0; 
        NV_DATA_S(ydot)[0] =  - ReactionFlux10 - ReactionFlux14 + ReactionFlux15 - ReactionFlux16 + ReactionFlux17;
    }    

    //d(V_C.T1)/dt
    NV_DATA_S(ydot)[1] =  - ReactionFlux26 - ReactionFlux32 + ReactionFlux33 - ReactionFlux34 + ReactionFlux35;

    //d(V_C.nivo)/dt
    NV_DATA_S(ydot)[2] = 1/PARAM(0)*( - ReactionFlux67 - ReactionFlux68 - ReactionFlux69 + ReactionFlux71 - ReactionFlux72);

    //d(V_C.durv)/dt
    NV_DATA_S(ydot)[3] = 1/PARAM(0)*( - ReactionFlux73 - ReactionFlux74 - ReactionFlux75 + ReactionFlux77 - ReactionFlux78);

    //d(V_C.ipi)/dt
    NV_DATA_S(ydot)[4] = 1/PARAM(0)*( - ReactionFlux79 - ReactionFlux80 - ReactionFlux81 + ReactionFlux83 - ReactionFlux84);

    //d(V_C.ENT)/dt
    NV_DATA_S(ydot)[5] = 1/PARAM(0)*(ReactionFlux146 + ReactionFlux147 - ReactionFlux148 - ReactionFlux149 - ReactionFlux150 - ReactionFlux151 + ReactionFlux153 - ReactionFlux154);

    //d(V_C.ENT_Buccal)/dt
    NV_DATA_S(ydot)[6] = 1/PARAM(0)*( - ReactionFlux146);

    //d(V_C.ENT_GI)/dt
    NV_DATA_S(ydot)[7] = 1/PARAM(0)*(ReactionFlux145 - ReactionFlux147);

    //d(V_C.Dose2)/dt
    NV_DATA_S(ydot)[8] = 1/PARAM(0)*( - ReactionFlux145);

    //d(V_P.T0)/dt
    NV_DATA_S(ydot)[9] =  - ReactionFlux11 + ReactionFlux14 - ReactionFlux15 - ReactionFlux134;

    //d(V_P.T1)/dt
    NV_DATA_S(ydot)[10] =  - ReactionFlux27 + ReactionFlux32 - ReactionFlux33;

    //d(V_P.nivo)/dt
    NV_DATA_S(ydot)[11] = 1/PARAM(1)*(ReactionFlux67);

    //d(V_P.durv)/dt
    NV_DATA_S(ydot)[12] = 1/PARAM(1)*(ReactionFlux73);

    //d(V_P.ipi)/dt
    NV_DATA_S(ydot)[13] = 1/PARAM(1)*(ReactionFlux79);

    //d(V_P.Treg_CTLA4)/dt
    NV_DATA_S(ydot)[14] =  - ReactionFlux131 - ReactionFlux132;

    //d(V_P.Treg_CTLA4_ipi)/dt
    NV_DATA_S(ydot)[15] = ReactionFlux131 - ReactionFlux132;

    //d(V_P.Treg_CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[16] = ReactionFlux132;

    //d(V_P.ENT)/dt
    NV_DATA_S(ydot)[17] = 1/PARAM(1)*(ReactionFlux149);

    //d(V_T.C_x)/dt
    NV_DATA_S(ydot)[18] =  (- ReactionFlux1 + ReactionFlux4 + ReactionFlux37);

    //d(V_T.T_exh)/dt
    NV_DATA_S(ydot)[19] =  (- ReactionFlux2 + ReactionFlux12 + ReactionFlux28 + ReactionFlux30 + ReactionFlux31);

    //d(V_T.C1)/dt
    NV_DATA_S(ydot)[20] = (ReactionFlux3 - ReactionFlux4 - ReactionFlux37);

    //d(V_T.T0)/dt

     if (SPVAR(21) / (PARAM(178) * AUX_VAR_V_T) < 1){
        realtype ReactionFlux16 = PARAM(27) * AUX_VAR_V_T * SPVAR(0) * (1.0 - SPVAR(21) / (PARAM(178) * AUX_VAR_V_T));
        realtype ReactionFlux144 = PARAM(172) * SPVAR(21) * (1.0 - SPVAR(21) / (PARAM(178) * AUX_VAR_V_T)) * SPVAR(35) / (PARAM(176) + SPVAR(35)) * (1.0 - SPVAR(36) / (SPVAR(36) + PARAM(181)));
        NV_DATA_S(ydot)[21] = ( - ReactionFlux12 + ReactionFlux16 - ReactionFlux133 + ReactionFlux144);
    }
    else{
        realtype ReactionFlux16 = 0; 
        realtype ReactionFlux144 = 0;
        NV_DATA_S(ydot)[21] = ( - ReactionFlux12 + ReactionFlux16 - ReactionFlux133 + ReactionFlux144);
    }   

    //d(V_T.T1)/dt
    NV_DATA_S(ydot)[22] = ( - ReactionFlux28 - ReactionFlux30 - ReactionFlux31 + ReactionFlux34);

    //d(V_T.APC)/dt
    NV_DATA_S(ydot)[23] = ReactionFlux38 - ReactionFlux40;

    //d(V_T.mAPC)/dt
    NV_DATA_S(ydot)[24] = ReactionFlux40 - ReactionFlux41 - ReactionFlux42;

    //d(V_T.c)/dt
    NV_DATA_S(ydot)[25] = 1/AUX_VAR_V_T*(ReactionFlux44 + ReactionFlux45);

    //d(V_T.nivo)/dt
    NV_DATA_S(ydot)[26] = 1/AUX_VAR_V_T*(ReactionFlux68 - ReactionFlux70);

    //d(V_T.durv)/dt
    NV_DATA_S(ydot)[27] = 1/AUX_VAR_V_T*(ReactionFlux74 - ReactionFlux76);

    //d(V_T.ipi)/dt
    NV_DATA_S(ydot)[28] = 1/AUX_VAR_V_T*(ReactionFlux80 - ReactionFlux82);

    //d(V_T.Treg_CTLA4)/dt
    NV_DATA_S(ydot)[29] =  - ReactionFlux129 - ReactionFlux130;

    //d(V_T.Treg_CTLA4_ipi)/dt
    NV_DATA_S(ydot)[30] = ReactionFlux129 - ReactionFlux130;

    //d(V_T.Treg_CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[31] = ReactionFlux130;

    //d(V_T.MDSC)/dt
    if (SPVAR(32) < PARAM(177) * AUX_VAR_V_T){

        realtype ReactionFlux135 = PARAM(162) * (PARAM(177) * AUX_VAR_V_T - SPVAR(32)) * (SPVAR(33) / (PARAM(175) + SPVAR(33)));

        realtype ReactionFlux136 = PARAM(180) * (PARAM(177) * AUX_VAR_V_T - SPVAR(32));

        NV_DATA_S(ydot)[32] = (ReactionFlux135 + ReactionFlux136 - ReactionFlux137);
    }
    else{

        realtype ReactionFlux135 = 0;

        realtype ReactionFlux136 = 0;

        NV_DATA_S(ydot)[32] = (ReactionFlux135 + ReactionFlux136 - ReactionFlux137);       
    }

    //d(V_T.CCL2)/dt
    NV_DATA_S(ydot)[33] = 1/AUX_VAR_V_T*( - ReactionFlux138 + ReactionFlux141);

    //d(V_T.NO)/dt
    NV_DATA_S(ydot)[34] = 1/AUX_VAR_V_T*( - ReactionFlux139 + ReactionFlux142);

    //d(V_T.ArgI)/dt
    NV_DATA_S(ydot)[35] = 1/AUX_VAR_V_T*( - ReactionFlux140 + ReactionFlux143);

    //d(V_T.ENT)/dt
    NV_DATA_S(ydot)[36] = 1/AUX_VAR_V_T*(ReactionFlux150 - ReactionFlux152);

    //d(V_LN.nT0)/dt
    NV_DATA_S(ydot)[37] = ReactionFlux5 - ReactionFlux6 - ReactionFlux7;

    //d(V_LN.aT0)/dt
    NV_DATA_S(ydot)[38] = ReactionFlux7 - ReactionFlux8;

    //d(V_LN.T0)/dt
    NV_DATA_S(ydot)[39] = ReactionFlux9 - ReactionFlux13 - ReactionFlux17;

    //d(V_LN.IL2)/dt
    NV_DATA_S(ydot)[40] = 1/PARAM(2)*( - ReactionFlux18 - ReactionFlux19 - ReactionFlux20 + ReactionFlux36);

    //d(V_LN.nT1)/dt
    NV_DATA_S(ydot)[41] = ReactionFlux21 - ReactionFlux22 - ReactionFlux23;

    //d(V_LN.aT1)/dt
    NV_DATA_S(ydot)[42] = ReactionFlux23 - ReactionFlux24;

    //d(V_LN.T1)/dt
    NV_DATA_S(ydot)[43] = ReactionFlux25 - ReactionFlux29 - ReactionFlux35;

    //d(V_LN.APC)/dt
    NV_DATA_S(ydot)[44] = ReactionFlux39;

    //d(V_LN.mAPC)/dt
    NV_DATA_S(ydot)[45] = ReactionFlux41 - ReactionFlux43;

    //d(V_LN.P0)/dt
    NV_DATA_S(ydot)[46] = 1/PARAM(2)*(ReactionFlux47 - ReactionFlux48 - ReactionFlux49);

    //d(V_LN.P1)/dt
    NV_DATA_S(ydot)[47] = 1/PARAM(2)*(ReactionFlux57 - ReactionFlux58 - ReactionFlux59);

    //d(V_LN.nivo)/dt
    NV_DATA_S(ydot)[48] = 1/PARAM(2)*(ReactionFlux69 + ReactionFlux70 - ReactionFlux71);

    //d(V_LN.durv)/dt
    NV_DATA_S(ydot)[49] = 1/PARAM(2)*(ReactionFlux75 + ReactionFlux76 - ReactionFlux77);

    //d(V_LN.ipi)/dt
    NV_DATA_S(ydot)[50] = 1/PARAM(2)*(ReactionFlux81 + ReactionFlux82 - ReactionFlux83);

    //d(V_LN.ENT)/dt
    NV_DATA_S(ydot)[51] = 1/PARAM(2)*(ReactionFlux151 + ReactionFlux152 - ReactionFlux153);

    //d(V_e.P0)/dt
    NV_DATA_S(ydot)[52] = 1/PARAM(3)*(ReactionFlux50 - ReactionFlux51);

    //d(V_e.p0)/dt
    NV_DATA_S(ydot)[53] = 1/PARAM(3)*(ReactionFlux51 - ReactionFlux52 - ReactionFlux53 + ReactionFlux54);

    //d(V_e.P1)/dt
    NV_DATA_S(ydot)[54] = 1/PARAM(3)*(ReactionFlux60 - ReactionFlux61);

    //d(V_e.p1)/dt
    NV_DATA_S(ydot)[55] = 1/PARAM(3)*(ReactionFlux61 - ReactionFlux62 - ReactionFlux63 + ReactionFlux64);

    //d(A_e.M1)/dt
    NV_DATA_S(ydot)[56] = 1/PARAM(4)*( - ReactionFlux46 - ReactionFlux53 + ReactionFlux54 - ReactionFlux63 + ReactionFlux64);

    //d(A_e.M1p0)/dt
    NV_DATA_S(ydot)[57] = 1/PARAM(4)*(ReactionFlux53 - ReactionFlux54 - ReactionFlux56);

    //d(A_e.M1p1)/dt
    NV_DATA_S(ydot)[58] = 1/PARAM(4)*(ReactionFlux63 - ReactionFlux64 - ReactionFlux66);

    //d(A_s.M1)/dt
    NV_DATA_S(ydot)[59] = 1/PARAM(5)*(ReactionFlux46 + ReactionFlux55 + ReactionFlux65);

    //d(A_s.M1p0)/dt
    NV_DATA_S(ydot)[60] = 1/PARAM(5)*( - ReactionFlux55 + ReactionFlux56);

    //d(A_s.M1p1)/dt
    NV_DATA_S(ydot)[61] = 1/PARAM(5)*( - ReactionFlux65 + ReactionFlux66);

    //d(syn_T_C1.PD1_PDL1)/dt
    NV_DATA_S(ydot)[62] = 1/PARAM(6)*(ReactionFlux85);

    //d(syn_T_C1.PD1_PDL2)/dt
    NV_DATA_S(ydot)[63] = 1/PARAM(6)*(ReactionFlux86);

    //d(syn_T_C1.PD1)/dt
    NV_DATA_S(ydot)[64] = 1/PARAM(6)*( - ReactionFlux85 - ReactionFlux86 - ReactionFlux87 - ReactionFlux88);

    //d(syn_T_C1.PDL1)/dt
    NV_DATA_S(ydot)[65] = 1/PARAM(6)*( - ReactionFlux85 - ReactionFlux89 - ReactionFlux90);

    //d(syn_T_C1.PDL2)/dt
    NV_DATA_S(ydot)[66] = 1/PARAM(6)*( - ReactionFlux86);

    //d(syn_T_C1.PD1_nivo)/dt
    NV_DATA_S(ydot)[67] = 1/PARAM(6)*(ReactionFlux87 - ReactionFlux88);

    //d(syn_T_C1.PD1_nivo_PD1)/dt
    NV_DATA_S(ydot)[68] = 1/PARAM(6)*(ReactionFlux88);

    //d(syn_T_C1.PDL1_durv)/dt
    NV_DATA_S(ydot)[69] = 1/PARAM(6)*(ReactionFlux89 - ReactionFlux90);

    //d(syn_T_C1.PDL1_durv_PDL1)/dt
    NV_DATA_S(ydot)[70] = 1/PARAM(6)*(ReactionFlux90);

    //d(syn_T_C1.TPDL1)/dt
    NV_DATA_S(ydot)[71] = 1/PARAM(6)*( - ReactionFlux103 - ReactionFlux104 - ReactionFlux105 - ReactionFlux106);

    //d(syn_T_C1.TPDL1_durv)/dt
    NV_DATA_S(ydot)[72] = 1/PARAM(6)*(ReactionFlux105 - ReactionFlux106);

    //d(syn_T_C1.TPDL1_durv_TPDL1)/dt
    NV_DATA_S(ydot)[73] = 1/PARAM(6)*(ReactionFlux106);

    //d(syn_T_C1.CD28_CD80)/dt
    NV_DATA_S(ydot)[74] = 1/PARAM(6)*(ReactionFlux91 - ReactionFlux92);

    //d(syn_T_C1.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[75] = 1/PARAM(6)*(ReactionFlux92);

    //d(syn_T_C1.CD28_CD86)/dt
    NV_DATA_S(ydot)[76] = 1/PARAM(6)*(ReactionFlux93);

    //d(syn_T_C1.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[77] = 1/PARAM(6)*(ReactionFlux94 - ReactionFlux95 - ReactionFlux97);

    //d(syn_T_C1.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[78] = 1/PARAM(6)*(ReactionFlux97 - ReactionFlux98);

    //d(syn_T_C1.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[79] = 1/PARAM(6)*(ReactionFlux95 - ReactionFlux96);

    //d(syn_T_C1.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[80] = 1/PARAM(6)*(ReactionFlux96 + ReactionFlux98);

    //d(syn_T_C1.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[81] = 1/PARAM(6)*(ReactionFlux99 - ReactionFlux100);

    //d(syn_T_C1.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[82] = 1/PARAM(6)*(ReactionFlux100);

    //d(syn_T_C1.TPDL1_CD80)/dt
    NV_DATA_S(ydot)[83] = 1/PARAM(6)*(ReactionFlux103 - ReactionFlux104);

    //d(syn_T_C1.TPDL1_CD80_TPDL1)/dt
    NV_DATA_S(ydot)[84] = 1/PARAM(6)*(ReactionFlux104);

    //d(syn_T_C1.CD28)/dt
    NV_DATA_S(ydot)[85] = 1/PARAM(6)*( - ReactionFlux91 - ReactionFlux92 - ReactionFlux93);

    //d(syn_T_C1.CTLA4)/dt
    NV_DATA_S(ydot)[86] = 1/PARAM(6)*( - ReactionFlux94 - ReactionFlux95 - ReactionFlux98 - ReactionFlux99 - ReactionFlux101 - ReactionFlux102);

    //d(syn_T_C1.CD80)/dt
    NV_DATA_S(ydot)[87] = 1/PARAM(6)*( - ReactionFlux91 - ReactionFlux94 - ReactionFlux96 - ReactionFlux97 - ReactionFlux103);

    //d(syn_T_C1.CD86)/dt
    NV_DATA_S(ydot)[88] = 1/PARAM(6)*( - ReactionFlux93 - ReactionFlux99 - ReactionFlux100);

    //d(syn_T_C1.CTLA4_ipi)/dt
    NV_DATA_S(ydot)[89] = 1/PARAM(6)*(ReactionFlux101 - ReactionFlux102);

    //d(syn_T_C1.CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[90] = 1/PARAM(6)*(ReactionFlux102);

    //d(syn_T_APC.PD1_PDL1)/dt
    NV_DATA_S(ydot)[91] = 1/PARAM(7)*(ReactionFlux107);

    //d(syn_T_APC.PD1_PDL2)/dt
    NV_DATA_S(ydot)[92] = 1/PARAM(7)*(ReactionFlux108);

    //d(syn_T_APC.PD1)/dt
    NV_DATA_S(ydot)[93] = 1/PARAM(7)*( - ReactionFlux107 - ReactionFlux108 - ReactionFlux109 - ReactionFlux110);

    //d(syn_T_APC.PDL1)/dt
    NV_DATA_S(ydot)[94] = 1/PARAM(7)*( - ReactionFlux107 - ReactionFlux111 - ReactionFlux112);

    //d(syn_T_APC.PDL2)/dt
    NV_DATA_S(ydot)[95] = 1/PARAM(7)*( - ReactionFlux108);

    //d(syn_T_APC.PD1_nivo)/dt
    NV_DATA_S(ydot)[96] = 1/PARAM(7)*(ReactionFlux109 - ReactionFlux110);

    //d(syn_T_APC.PD1_nivo_PD1)/dt
    NV_DATA_S(ydot)[97] = 1/PARAM(7)*(ReactionFlux110);

    //d(syn_T_APC.PDL1_durv)/dt
    NV_DATA_S(ydot)[98] = 1/PARAM(7)*(ReactionFlux111 - ReactionFlux112);

    //d(syn_T_APC.PDL1_durv_PDL1)/dt
    NV_DATA_S(ydot)[99] = 1/PARAM(7)*(ReactionFlux112);

    //d(syn_T_APC.TPDL1)/dt
    NV_DATA_S(ydot)[100] = 1/PARAM(7)*( - ReactionFlux125 - ReactionFlux126 - ReactionFlux127 - ReactionFlux128);

    //d(syn_T_APC.TPDL1_durv)/dt
    NV_DATA_S(ydot)[101] = 1/PARAM(7)*(ReactionFlux127 - ReactionFlux128);

    //d(syn_T_APC.TPDL1_durv_TPDL1)/dt
    NV_DATA_S(ydot)[102] = 1/PARAM(7)*(ReactionFlux128);

    //d(syn_T_APC.CD28_CD80)/dt
    NV_DATA_S(ydot)[103] = 1/PARAM(7)*(ReactionFlux113 - ReactionFlux114);

    //d(syn_T_APC.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[104] = 1/PARAM(7)*(ReactionFlux114);

    //d(syn_T_APC.CD28_CD86)/dt
    NV_DATA_S(ydot)[105] = 1/PARAM(7)*(ReactionFlux115);

    //d(syn_T_APC.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[106] = 1/PARAM(7)*(ReactionFlux116 - ReactionFlux117 - ReactionFlux119);

    //d(syn_T_APC.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[107] = 1/PARAM(7)*(ReactionFlux119 - ReactionFlux120);

    //d(syn_T_APC.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[108] = 1/PARAM(7)*(ReactionFlux117 - ReactionFlux118);

    //d(syn_T_APC.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[109] = 1/PARAM(7)*(ReactionFlux118 + ReactionFlux120);

    //d(syn_T_APC.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[110] = 1/PARAM(7)*(ReactionFlux121 - ReactionFlux122);

    //d(syn_T_APC.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[111] = 1/PARAM(7)*(ReactionFlux122);

    //d(syn_T_APC.TPDL1_CD80)/dt
    NV_DATA_S(ydot)[112] = 1/PARAM(7)*(ReactionFlux125 - ReactionFlux126);

    //d(syn_T_APC.TPDL1_CD80_TPDL1)/dt
    NV_DATA_S(ydot)[113] = 1/PARAM(7)*(ReactionFlux126);

    //d(syn_T_APC.CD28)/dt
    NV_DATA_S(ydot)[114] = 1/PARAM(7)*( - ReactionFlux113 - ReactionFlux114 - ReactionFlux115);

    //d(syn_T_APC.CTLA4)/dt
    NV_DATA_S(ydot)[115] = 1/PARAM(7)*( - ReactionFlux116 - ReactionFlux117 - ReactionFlux120 - ReactionFlux121 - ReactionFlux123 - ReactionFlux124);

    //d(syn_T_APC.CD80)/dt
    NV_DATA_S(ydot)[116] = 1/PARAM(7)*( - ReactionFlux113 - ReactionFlux116 - ReactionFlux118 - ReactionFlux119 - ReactionFlux125);

    //d(syn_T_APC.CD86)/dt
    NV_DATA_S(ydot)[117] = 1/PARAM(7)*( - ReactionFlux115 - ReactionFlux121 - ReactionFlux122);

    //d(syn_T_APC.CTLA4_ipi)/dt
    NV_DATA_S(ydot)[118] = 1/PARAM(7)*(ReactionFlux123 - ReactionFlux124);

    //d(syn_T_APC.CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[119] = 1/PARAM(7)*(ReactionFlux124);

    return(0);
}
int ODE_system::g(realtype t, N_Vector y, realtype *gout, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    //V_T.C1 < (0.9 * cell)
    gout[0] = 0.9 * PARAM(9) - (SPVAR(20));

    return(0);
}

bool ODE_system::triggerComponentEvaluate(int i, realtype t, bool curr) {

    bool discrete = false;
    realtype diff = 0;
    bool eval = false;
    //Assignment rules:

    switch(i)
    {
    case 0:
        //V_T.C1 < (0.9 * cell)
        diff = 0.9 * _class_parameter[9] - (NV_DATA_S(_y)[20]);
        break;
    default:
        break;
    }
    if (!discrete){
        eval = diff == 0 ? curr : (diff > 0);
    }
    return eval;
}

bool ODE_system::eventEvaluate(int i) {
    bool eval = false;
    switch(i)
    {
    case 0:
        eval = getSatisfied(0);
        break;
    default:
        break;
    }
    return eval;
}

bool ODE_system::eventExecution(int i, bool delayed, realtype& dt){

    bool setDelay = false;

    //Assignment rules:

    switch(i)
    {
    case 0:
        NV_DATA_S(_y)[20] = 0.0 * _class_parameter[9];
        break;
    default:
        break;
    }
    return setDelay;
}
void ODE_system::update_y_other(void){

    return;
}
std::string ODE_system::getHeader(){

    std::string s = "";
    s += ",V_C.T0";
    s += ",V_C.T1";
    s += ",V_C.nivo";
    s += ",V_C.durv";
    s += ",V_C.ipi";
    s += ",V_C.ENT";
    s += ",V_C.ENT_Buccal";
    s += ",V_C.ENT_GI";
    s += ",V_C.Dose2";
    s += ",V_P.T0";
    s += ",V_P.T1";
    s += ",V_P.nivo";
    s += ",V_P.durv";
    s += ",V_P.ipi";
    s += ",V_P.Treg_CTLA4";
    s += ",V_P.Treg_CTLA4_ipi";
    s += ",V_P.Treg_CTLA4_ipi_CTLA4";
    s += ",V_P.ENT";
    s += ",V_T.C_x";
    s += ",V_T.T_exh";
    s += ",V_T.C1";
    s += ",V_T.T0";
    s += ",V_T.T1";
    s += ",V_T.APC";
    s += ",V_T.mAPC";
    s += ",V_T.c";
    s += ",V_T.nivo";
    s += ",V_T.durv";
    s += ",V_T.ipi";
    s += ",V_T.Treg_CTLA4";
    s += ",V_T.Treg_CTLA4_ipi";
    s += ",V_T.Treg_CTLA4_ipi_CTLA4";
    s += ",V_T.MDSC";
    s += ",V_T.CCL2";
    s += ",V_T.NO";
    s += ",V_T.ArgI";
    s += ",V_T.ENT";
    s += ",V_LN.nT0";
    s += ",V_LN.aT0";
    s += ",V_LN.T0";
    s += ",V_LN.IL2";
    s += ",V_LN.nT1";
    s += ",V_LN.aT1";
    s += ",V_LN.T1";
    s += ",V_LN.APC";
    s += ",V_LN.mAPC";
    s += ",V_LN.P0";
    s += ",V_LN.P1";
    s += ",V_LN.nivo";
    s += ",V_LN.durv";
    s += ",V_LN.ipi";
    s += ",V_LN.ENT";
    s += ",V_e.P0";
    s += ",V_e.p0";
    s += ",V_e.P1";
    s += ",V_e.p1";
    s += ",A_e.M1";
    s += ",A_e.M1p0";
    s += ",A_e.M1p1";
    s += ",A_s.M1";
    s += ",A_s.M1p0";
    s += ",A_s.M1p1";
    s += ",syn_T_C1.PD1_PDL1";
    s += ",syn_T_C1.PD1_PDL2";
    s += ",syn_T_C1.PD1";
    s += ",syn_T_C1.PDL1";
    s += ",syn_T_C1.PDL2";
    s += ",syn_T_C1.PD1_nivo";
    s += ",syn_T_C1.PD1_nivo_PD1";
    s += ",syn_T_C1.PDL1_durv";
    s += ",syn_T_C1.PDL1_durv_PDL1";
    s += ",syn_T_C1.TPDL1";
    s += ",syn_T_C1.TPDL1_durv";
    s += ",syn_T_C1.TPDL1_durv_TPDL1";
    s += ",syn_T_C1.CD28_CD80";
    s += ",syn_T_C1.CD28_CD80_CD28";
    s += ",syn_T_C1.CD28_CD86";
    s += ",syn_T_C1.CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80";
    s += ",syn_T_C1.CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4_CD86";
    s += ",syn_T_C1.TPDL1_CD80";
    s += ",syn_T_C1.TPDL1_CD80_TPDL1";
    s += ",syn_T_C1.CD28";
    s += ",syn_T_C1.CTLA4";
    s += ",syn_T_C1.CD80";
    s += ",syn_T_C1.CD86";
    s += ",syn_T_C1.CTLA4_ipi";
    s += ",syn_T_C1.CTLA4_ipi_CTLA4";
    s += ",syn_T_APC.PD1_PDL1";
    s += ",syn_T_APC.PD1_PDL2";
    s += ",syn_T_APC.PD1";
    s += ",syn_T_APC.PDL1";
    s += ",syn_T_APC.PDL2";
    s += ",syn_T_APC.PD1_nivo";
    s += ",syn_T_APC.PD1_nivo_PD1";
    s += ",syn_T_APC.PDL1_durv";
    s += ",syn_T_APC.PDL1_durv_PDL1";
    s += ",syn_T_APC.TPDL1";
    s += ",syn_T_APC.TPDL1_durv";
    s += ",syn_T_APC.TPDL1_durv_TPDL1";
    s += ",syn_T_APC.CD28_CD80";
    s += ",syn_T_APC.CD28_CD80_CD28";
    s += ",syn_T_APC.CD28_CD86";
    s += ",syn_T_APC.CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80";
    s += ",syn_T_APC.CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4_CD86";
    s += ",syn_T_APC.TPDL1_CD80";
    s += ",syn_T_APC.TPDL1_CD80_TPDL1";
    s += ",syn_T_APC.CD28";
    s += ",syn_T_APC.CTLA4";
    s += ",syn_T_APC.CD80";
    s += ",syn_T_APC.CD86";
    s += ",syn_T_APC.CTLA4_ipi";
    s += ",syn_T_APC.CTLA4_ipi_CTLA4";
    return s;
}
realtype ODE_system::get_unit_conversion_species(int i) const{

    static std::vector<realtype> scalor = {
        //sp_var
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000000.0,
        1000000.0,
        1000000.0,
        1000000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000000.0,
        1000000.0,
        1000000.0,
        1000000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999999999.9999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        //sp_other
    };
    return scalor[i];
}
realtype ODE_system::get_unit_conversion_nspvar(int i) const{

    static std::vector<realtype> scalor = {
    };
    return scalor[i];
}
};
