using Random
using CartesianGeneticProgramming
include("../scripts/model_template_utils.jl")

function get_test_ind(cfg::NamedTuple)
    local inputs, template, outputs
    inputs = [
        "x1",
        "x2",
        "x3",
    ]
    template = [
        (id="f1", f="f_mult", x="x1", y="x2"),
        (id="self+f1", f="f_add", x="self+f1", y="f1")
    ]
    outputs = ["self+f1"]
    return get_ind_from_template(cfg, template, inputs, outputs)
end

function get_shorefor(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_P",#2
        "phi",#3
        "0.5#1",#4 
        "0.5#2",#5
        "dt_waves",#6 
        "current_time",#7 
        "rand",#8 
        "coeff_time",#9
        "coeff_v",#10
        "r"#11
    ]

    shorefor = [
        (id="indices", f="f_irange", x="input_omega", y=""),#12
        (id="flipped_indices", f="f_reverse", x="indices", y=""),#22
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),#32
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),#42
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),#52
        (id="weighted_omega", f="f_mult", x="input_omega", y="weights"),#62
        (id="sum_weights", f="f_sum", x="weights", y=""),#72
        (id="sum_weighted_omega", f="f_sum", x="weighted_omega", y=""),#82
        (id="omega_eq", f="f_div", x="sum_weighted_omega", y="sum_weights"),#92
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),#102
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),#112
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),#122
        (id="input_P^0.5", f="f_pow", x="input_P", y="0.5#1", p=["0.5#1"]),#132
        (id="F", f="f_mult", x="input_P^0.5", y="norm_d_omega"),#142
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),#152
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),#162
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),#172
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),#182
        (id="s1", f="f_mult", x="Fe", y="r", p=["r"]),#192
        (id="s2", f="f_add", x="s1", y="Fa"),#202
        (id="ms_s2", f="f_moving_sum", x="s2", y=""),#212
        (id="end_ms_s2", f="f_pop", x="ms_s2", y=""),#222
        (id="half_end_ms_s2", f="f_mult", x="end_ms_s2", y="0.5#2", p=["0.5#2"]),#232
        (id="dtw_half_end_ms_s2", f="f_mult", x="half_end_ms_s2", y="dt_waves"),#242
        (id="v", f="f_add", x="v", y="dtw_half_end_ms_s2"),#252
        (id="sx1", f="f_mult", x="coeff_time", y="current_time"),#262
        (id="sx2", f="f_mult", x="coeff_v", y="v"),#272
        # (id="sx3", f="f_add", x="sx1", y="sx2"),#282
        # (id="x", f="f_add", x="sx3", y="intercept", p=["intercept"]),#292
        (id="x", f="f_add", x="sx1", y="sx2")#282
    ]

    outputs = ["x"]

    # shorefor_ind = get_ind_from_template(cfg, shorefor, inputs, outputs)
    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_shorefor_dxdt(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_P",#2
        "phi",#3
        "0.5#1",#4 
        "r",#5
        "input_dir"
    ]

    shorefor = [
        (id="indices", f="f_irange", x="input_omega", y=""),#12
        (id="flipped_indices", f="f_reverse", x="indices", y=""),#22
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),#32
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),#42
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),#52
        (id="weighted_omega", f="f_mult", x="input_omega", y="weights"),#62
        (id="sum_weights", f="f_sum", x="weights", y=""),#72
        (id="sum_weighted_omega", f="f_sum", x="weighted_omega", y=""),#82
        (id="omega_eq", f="f_div", x="sum_weighted_omega", y="sum_weights"),#92
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),#102
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),#112
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),#122
        (id="input_P^0.5", f="f_pow", x="input_P", y="0.5#1", p=["0.5#1"]),#132
        (id="F", f="f_mult", x="input_P^0.5", y="norm_d_omega"),#142
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),#152
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),#162
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),#172
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),#182
        (id="s1", f="f_mult", x="Fe", y="r", p=["r"]),#192
        (id="s2", f="f_add", x="s1", y="Fa"),#202
        (id="dxdt", f="f_pop", x="s2", y="")#282
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_fullmodel_shorefor_dxdt(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_omega_eq",#2
        "input_P",#3
        "phi",#4
        "0.5#1",#5
        "Dir"#6
    ]

    shorefor = [
        (id="d_omega", f="f_subtract", x="input_omega_eq", y="input_omega"),
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),
        (id="input_P^0.5", f="f_pow", x="input_P", y="0.5#1", p=["0.5#1"]),
        (id="F", f="f_mult", x="input_P^0.5", y="norm_d_omega"),
        (id="erosion_bools", f="f_lteq", x="input_omega_eq", y="input_omega"),
        (id="accretion_bools", f="f_gt", x="input_omega_eq", y="input_omega"),
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),
        (id="integrate_Fa", f="f_integrate", x="Fa", y=""),
        (id="integrate_Fe", f="f_integrate", x="Fe", y=""),
        (id="Fa/Fe", f="f_div", x="integrate_Fa", y="integrate_Fe"),
        (id="calc_r", f="f_abs", x="Fa/Fe", y=""),
        (id="s1", f="f_mult", x="Fe", y="calc_r"),
        (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["dxdt"]#, "calc_r", "Fa", "Fe"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_fullmodel_shorefor_dxdt_v2(cfg::NamedTuple; return_ids_indices=false)
    local shorefor, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]

    shorefor = [
        (id="indices", f="f_irange", x="2phi", y=""),
        (id="flipped_indices", f="f_reverse", x="indices", y=""),
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),
        (id="sum_weights", f="f_sum", x="weights", y=""),
        (id="weights_filter", f="f_div", x="weights", y="sum_weights"),
        (id="omega_eq", f="f_conv", x="input_omega", y="weights_filter"),
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),
        (id="F", f="f_mult", x="input_{P^0.5}", y="norm_d_omega"),
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),
        (id="detrend_Fa", f="f_detrend", x="Fa", y=""),
        (id="detrend_Fe", f="f_detrend", x="Fe", y=""),
        (id="integrate_Fa", f="f_integrate", x="detrend_Fa", y=""),
        (id="integrate_Fe", f="f_integrate", x="detrend_Fe", y=""),
        (id="Fa/Fe", f="f_div", x="integrate_Fa", y="integrate_Fe"),
        (id="calc_r", f="f_abs", x="Fa/Fe", y=""),
        (id="s1", f="f_mult", x="Fe", y="calc_r"),
        (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_recurrent_shorefor_dxdt(cfg::NamedTuple; return_ids_indices=false)
    local shorefor, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]

    shorefor = [
        (id="indices", f="f_irange", x="2phi", y=""),
        (id="flipped_indices", f="f_reverse", x="indices", y=""),
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),
        (id="sum_weights", f="f_sum", x="weights", y=""),
        (id="weights_filter", f="f_div", x="weights", y="sum_weights"),
        (id="omega_eq", f="f_conv", x="input_omega", y="weights_filter"),
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),
        (id="F", f="f_mult", x="input_{P^0.5}", y="norm_d_omega"),
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),
        (id="detrend_Fa", f="f_detrend", x="Fa", y=""),
        (id="detrend_Fe", f="f_detrend", x="Fe", y=""),
        (id="integrate_Fa", f="f_integrate", x="detrend_Fa", y=""),
        (id="integrate_Fe", f="f_integrate", x="detrend_Fe", y=""),
        (id="Fa/Fe", f="f_div", x="integrate_Fa", y="integrate_Fe"),
        (id="calc_r", f="f_abs", x="Fa/Fe", y=""),
        (id="s1", f="f_mult", x="Fe", y="calc_r"),
        (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_shorefor_evolved(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_P",#2
        "phi",#3
        "0.5#1",#4 
        "0.5#2",#5
        "dt_waves",#6 
        "current_time",#7 
        "rand",#8 
        "coeff_time",#9
        "coeff_v",#10
        "r"#11
    ]

    shorefor = [
        (id="indices", f="f_irange", x="input_omega", y=""),#12
        (id="flipped_indices", f="f_pop", x="indices", y=""),#22  #f_pop    ******
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),#32
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="indices", p=[""]),#42  #y=indices   ******
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),#52
        (id="weighted_omega", f="f_mult", x="input_omega", y="weights"),#62
        (id="sum_weights", f="f_sum", x="weights", y=""),#72
        (id="sum_weighted_omega", f="f_sum", x="weighted_omega", y=""),#82
        (id="omega_eq", f="f_div", x="sum_weighted_omega", y="sum_weights"),#92
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),#102
        # (id="std_d_omega", f="f_stddev", x="d_omega", y=""),#112  ## REMOVED NORM TERM
        (id="norm_d_omega", f="f_moving_sum", x="d_omega", y=""),#112
        (id="input_P^0.5", f="f_pow", x="input_P", y="0.5#1", p=["0.5#1"]),#122
        (id="F", f="f_mult", x="input_P^0.5", y="norm_d_omega"),#132
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),#142
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),#152
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),#162
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),#172
        (id="s1", f="f_mult", x="Fe", y="r", p=["r"]),#182
        (id="s2", f="f_add", x="s1", y="Fa"),#192
        (id="ms_s2", f="f_moving_sum", x="s2", y=""),#202
        (id="end_ms_s2", f="f_pop", x="ms_s2", y=""),#212
        (id="half_end_ms_s2", f="f_mult", x="end_ms_s2", y="coeff_time", p=["0.5#2"]),#222
        (id="dtw_half_end_ms_s2", f="f_mult", x="half_end_ms_s2", y="dt_waves"),#232
        (id="v", f="f_add", x="v", y="dtw_half_end_ms_s2"),#242
    ]

    outputs = ["v"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_module_omegaeq_normdeltaomega(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",
        "phi",
    ]

    shorefor = [
        (id="indices", f="f_irange", x="input_omega", y=""),#12
        (id="flipped_indices", f="f_reverse", x="indices", y=""),#22
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),#32
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),#42
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),#52
        (id="weighted_omega", f="f_mult", x="input_omega", y="weights"),#62
        (id="sum_weights", f="f_sum", x="weights", y=""),#72
        (id="sum_weighted_omega", f="f_sum", x="weighted_omega", y=""),#82
        (id="omega_eq", f="f_div", x="sum_weighted_omega", y="sum_weights"),#92
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),#102
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),#112
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),#122
    ]

    outputs = ["omega_eq", "norm_d_omega"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_module_omegaeq(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",
        "phi",
        "input_Dir"
    ]

    shorefor = [
        (id="indices", f="f_irange", x="input_omega", y=""),#12
        (id="flipped_indices", f="f_reverse", x="indices", y=""),#22
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),#32
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),#42
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),#52
        (id="weighted_omega", f="f_mult", x="input_omega", y="weights"),#62
        (id="sum_weights", f="f_sum", x="weights", y=""),#72
        (id="sum_weighted_omega", f="f_sum", x="weighted_omega", y=""),#82
        (id="omega_eq", f="f_div", x="sum_weighted_omega", y="sum_weights"),
        (id="dir_eq", f="f_mean", x="input_Dir", y="")
    ]

    outputs = ["omega_eq", "dir_eq"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_module_dxdt(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_omega_eq",#2
        "input_P",#3
        "phi",#4
        "0.5#1",#5
        "input_Dir",#6,
        "input_Dir_eq"#7
    ]

    shorefor = [
        (id="d_omega", f="f_subtract", x="input_omega_eq", y="input_omega"),
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),
        (id="input_P^0.5", f="f_pow", x="input_P", y="0.5#1", p=["0.5#1"]),
        (id="F", f="f_mult", x="input_P^0.5", y="norm_d_omega"),
        (id="erosion_bools", f="f_lteq", x="input_omega_eq", y="input_omega"),
        (id="accretion_bools", f="f_gt", x="input_omega_eq", y="input_omega"),
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),
        (id="integrate_Fa", f="f_integrate", x="Fa", y=""),
        (id="integrate_Fe", f="f_integrate", x="Fe", y=""),
        (id="Fa/Fe", f="f_div", x="integrate_Fa", y="integrate_Fe"),
        (id="calc_r", f="f_abs", x="Fa/Fe", y=""),
        (id="s1", f="f_mult", x="Fe", y="calc_r"),
        (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_shorefor_evolved_fmconv_oldins(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
    ]

    # shorefor = [
    #     (id="tail_tp", f="f_tail", x="Tp", y=""),
    #     (id="vec_from_tail_tp", f="f_vecfromdouble", x="tail_tp", y=""),
    #     (id="sin_hsb", f="f_sin", x="Hsb", y=""),
    #     (id="diff_sin_hsb", f="f_diff", x="sin_hsb", y=""),
    #     (id="conv1", f="f_conv", x="vec_from_tail_tp", y="diff_sin_hsb"),
    #     (id="ceil_conv1", f="f_ceil", x="conv1", y=""),
    #     (id="cos_p", f="f_cos", x="input_{P^0.5}", y=""),
    #     (id="norm_p", f="f_normalize", x="input_{P^0.5}", y=""),
    #     (id="lteq1", f="f_lteq", x="cos_p", y="ceil_conv1"),
    #     (id="min1", f="f_min", x="lteq1", y="norm_p"),
    #     (id="sqrt_min1", f="f_sqrt", x="min1", y=""),
    #     (id="ceil_recur", f="f_ceil", x="ceil_recur", y=""),
    #     (id="sqrtxy_dir_ceilrecur", f="f_sqrt_xy", x="Dir", y="ceil_recur"),
    #     (id="conv2", f="f_conv", x="phi", y="conv2"),
    #     (id="conv3", f="f_conv", x="sqrtxy_dir_ceilrecur", y="conv2"),
        
    #     (id="dxdt", f="f_mult", x="sqrt_min1", y="conv3"),
    #     ]

    shorefor = [
        (id="tail_tp", f="f_tail", x="Tp", y=""),
        # (id="vec_from_tail_tp", f="f_vecfromdouble", x="tail_tp", y=""),
        (id="sin_hsb", f="f_sin", x="Hsb", y=""),
        (id="diff_sin_hsb", f="f_diff", x="sin_hsb", y=""),
        (id="conv1", f="f_conv", x="tail_tp", y="diff_sin_hsb"),
        (id="ceil_conv1", f="f_ceil", x="conv1", y=""),
        (id="cos_p", f="f_cos", x="input_{P^0.5}", y=""),
        (id="norm_p", f="f_normalize", x="input_{P^0.5}", y=""),
        (id="lteq1", f="f_lteq", x="cos_p", y="ceil_conv1"),
        (id="min1", f="f_min", x="lteq1", y="norm_p"),
        (id="sqrt_min1", f="f_sqrt", x="min1", y=""),
        (id="sqrt_dir", f="f_sqrt", x="Dir", y=""),
        (id="dxdt", f="f_mult", x="sqrt_min1", y="sqrt_dir"),
        ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_vis_omegaeq(cfg::NamedTuple)
    local omega_eq, inputs, outputs

    inputs = [
        "input_omega",#1
        "phi",#3
        "2phi",#4
    ]

    omega_eq = [
        (id="indices", f="f_irange", x="2phi", y=""),
        (id="flipped_indices", f="f_reverse", x="indices", y=""),
        (id="neg_flip_indices", f="f_negate", x="flipped_indices", y=""),
        (id="neg_flip_indices/phi", f="f_div", x="neg_flip_indices", y="phi", p=["phi"]),
        (id="weights", f="f_tpow", x="neg_flip_indices/phi", y=""),
        (id="sum_weights", f="f_sum", x="weights", y=""),
        (id="weights_filter", f="f_div", x="weights", y="sum_weights"),
        (id="omega_eq", f="f_conv", x="input_omega", y="weights_filter"),
    ]

    outputs = ["omega_eq"]

    return get_ind_from_template(cfg, omega_eq, inputs, outputs)
end

function get_vis_F(cfg::NamedTuple)
    local shorefor, inputs, outputs

    inputs = [
        "input_omega",#1
        "omega_eq",
        "input_{P^0.5}"
    ]

    shorefor = [
        (id="d_omega", f="f_subtract", x="omega_eq", y="input_omega"),
        (id="std_d_omega", f="f_stddev", x="d_omega", y=""),
        (id="norm_d_omega", f="f_div", x="d_omega", y="std_d_omega"),
        (id="F", f="f_mult", x="input_{P^0.5}", y="norm_d_omega"),
    ]

    outputs = ["F"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_vis_r(cfg::NamedTuple; return_ids_indices=false)
    local shorefor, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "Fe",
        "Fa"
    ]

    shorefor = [
        (id="detrend_Fa", f="f_detrend", x="Fa", y=""),
        (id="detrend_Fe", f="f_detrend", x="Fe", y=""),
        (id="integrate_Fa", f="f_integrate", x="detrend_Fa", y=""),
        (id="integrate_Fe", f="f_integrate", x="detrend_Fe", y=""),
        (id="Fa/Fe", f="f_div", x="integrate_Fa", y="integrate_Fe"),
        (id="r", f="f_abs", x="Fa/Fe", y=""),
        # (id="s1", f="f_mult", x="Fe", y="calc_r"),
        # (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["r"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_vis_FeFa(cfg::NamedTuple; return_ids_indices=false)
    local shorefor, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "omega_eq",
        "F"
    ]

    shorefor = [
        
        (id="erosion_bools", f="f_lteq", x="omega_eq", y="input_omega"),
        (id="accretion_bools", f="f_gt", x="omega_eq", y="input_omega"),
        (id="Fe", f="f_mult", x="F", y="erosion_bools"),
        (id="Fa", f="f_mult", x="F", y="accretion_bools"),
    ]

    outputs = ["Fe", "Fa"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_vis_dxdt(cfg::NamedTuple; return_ids_indices=false)
    local shorefor, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "Fe",#1
        "Fa",#2
        "r",
    ]

    shorefor = [
        (id="s1", f="f_mult", x="Fe", y="r"),
        (id="dxdt", f="f_add", x="s1", y="Fa"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, shorefor, inputs, outputs)
end

function get_global_P_ind(cfg::NamedTuple; return_ids_indices=false)
    local model, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]

    model = [
        (id="10", f="f_pow", x="Tp", y="Dir"),
        (id="12", f="f_conv", x="10", y="Dir"),
        (id="13", f="f_sqrt", x="12", y=""),
        (id="16", f="f_sqrt_xy", x="Dir", y="13"),
        (id="18", f="f_sqrt", x="16", y=""),
        (id="22", f="f_sqrt_xy", x="16", y="18"),
        (id="21", f="f_add", x="Dir", y="input_omega"),
        (id="23", f="f_div", x="22", y="21"),
        (id="11", f="f_diff", x="10", y=""),
        (id="20", f="f_add", x="13", y="11"),
        (id="42", f="f_div", x="20", y="23"),
        (id="dxdt", f="f_reverse", x="42", y=""),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, model, inputs, outputs)
end

function get_modified_global_P_ind(cfg::NamedTuple; return_ids_indices=false)
    local model, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]

    model = [
        (id="10", f="f_pow", x="Tp", y="Dir"),
        (id="12", f="f_mean", x="10", y=""),
        (id="13", f="f_sqrt", x="12", y=""),
        (id="16", f="f_sqrt_xy", x="Dir", y="13"),
        (id="18", f="f_sqrt", x="16", y=""),
        (id="22", f="f_sqrt_xy", x="16", y="18"),
        (id="21", f="f_add", x="Dir", y="input_omega"),
        (id="23", f="f_div", x="22", y="21"),
        # (id="23", f="f_div", x="22", y="16"),
        (id="11", f="f_diff", x="10", y=""),
        (id="20", f="f_add", x="13", y="11"),
        # (id="dxdt", f="f_div", x="20", y="23"),
        (id="42", f="f_div", x="20", y="23"),
        (id="dxdt", f="f_reverse", x="42", y=""),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, model, inputs, outputs)
end

function get_global_P_4_ind(cfg::NamedTuple; return_ids_indices=false)
    local model, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]

    # inputs = [
    #     "input_omega",
    #     "Dir",
    #     "Hsb"
    # ]

    model = [
        # (id="24", f="f_mult", x="Dir", y="Hsb"),
        # (id="36", f="f_gt", x="24", y="phi"),
        # (id="38", f="f_nop", x="36", y=""),
        # (id="23", f="f_sqrt", x="input_omega", y=""),
        # (id="35", f="f_pow", x="24", y="23"),
        # (id="49", f="f_subtract", x="35", y="38"),
        # (id="dxdt", f="f_nop", x="49", y="")
        (id="24", f="f_mult", x="Dir", y="Hsb"),
        # (id="36", f="f_gt", x="24", y="phi"),
        (id="23", f="f_sqrt", x="input_omega", y=""),
        (id="dxdt", f="f_pow", x="24", y="23"),
        # (id="dxdt", f="f_subtract", x="35", y="36"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, model, inputs, outputs)
end

function get_global_SLA_ind(cfg::NamedTuple; return_ids_indices=false)
    local model, inputs, outputs

    # omega, P, phi, 2phi, Dir, Hsb, Tp, E, Sla, rivdis
    inputs = [
        "input_omega",#1
        "input_{P^0.5}",#2
        "phi",#3
        "2phi",#4
        "Dir",#5
        "Hsb",#6
        "Tp",#7
        "Sla",#9
        "rivdis"#10
    ]
    
    ## for chromo_draw:
    # inputs = [
    #     "input_omega",#1
    #     "Sla",#9
    # ]

    model = [
        (id="10", f="f_sin", x="Sla", y=""),
        (id="12", f="f_reverse", x="10", y=""),
        (id="14", f="f_add", x="12", y="10"),
        (id="17", f="f_pow", x="14", y="input_omega"),
        (id="18", f="f_tpow", x="17", y=""),
        (id="19", f="f_diff", x="18", y=""),
        (id="21", f="f_add", x="19", y="Sla"),
        (id="24", f="f_diff", x="21", y=""),
        (id="22", f="f_negate", x="19", y=""),
        (id="15", f="f_subtract", x="10", y="input_omega"),
        (id="16", f="f_subtract", x="15", y="Sla"),
        (id="29", f="f_subtract", x="16", y="22"),
        (id="32", f="f_add", x="29", y="24"),
        (id="49", f="f_reverse", x="32", y=""),
        (id="dxdt", f="f_add", x="49", y="24"),
    ]

    outputs = ["dxdt"]

    return get_ind_from_template(cfg, model, inputs, outputs)
end