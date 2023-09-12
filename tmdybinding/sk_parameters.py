from .tmd_abstract_lattice import SKSimpleParametersList, SKParametersList


_params_sk_cappelluti = SKSimpleParametersList(dict(zip(
    # from https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.075409
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.512,         None,           -3.025,         -1.276,         -8.236,
        -2.619,         -1.396,         -0.933,         -0.478,         -0.442,
        0.696,          0.278,          0.696,           0.278,
    ]
)))

_params_sk_roldan_mos2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/2053-1583/1/3/034003
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.512,         0.419,          -3.025,         -1.276,         -8.236,
        -2.619,         -1.396,         -0.933,         -0.478,         -0.442,
        0.696,          0.278,          0.696,          0.278
    ]
)))

_params_sk_roldan_ws2 = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/2053-1583/1/3/034003
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3153,         0.215,          0.057,          "WS2",          0.712,
        -1.550,         0.851,          -3.090,         -1.176,         -7.836,
        -0.619,         -1.396,         -0.983,         -0.478,         -0.442,
        0.696,          0.278,          0.696,          0.278
    ]
)))

_params_sk_ridolfi = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        0.201,          -1.563,         -0.352,         -54.839,        -39.275,
        -9.88,          4.196,          -1.153,         0.612,          0.086,
        12.734,         -2.175,         12.734,         -2.175
    ]
)))


_params_sk_ridolfi_vb = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        0.191,          -1.599,         0.081,          -48.934,        -37.981,
        -8.963,         4.115,          -1.154,         0.964,          0.117,
        10.707,         -4.084,         10.707,         -4.084
    ]
)))

_params_sk_ridolfi_minimal = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1088/0953-8984/27/36/365501
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.707,
        -11.683,        -208.435,       -75.952,        -23.761,        -35.968,
        -56.738,        1.318,          -2.652,         1.750,          1.482,
        0.000,          0.000,          0.000,          0.000
    ]
)))

_params_sk_venkateswarlu = SKSimpleParametersList(dict(zip(
    # from https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.081103
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.318,          0.086,          0.052,          "MoS2",         0.704,
        0.1356,         -0.4204,        0.0149,         -38.71,         -29.45,
        -7.193,         3.267,          -0.9035,         0.7027,        0.0897,
        8.079,          -2.678,         7.336,          -2.432
    ]
)))

_params_sk_silva_guillen_mos2 = SKSimpleParametersList(dict(zip(
    # from https://www.mdpi.com/2076-3417/6/10/284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.086,          0.052,          "MoS2",         0.716,
        -1.094,         -0.050,         -1.511,         -3.559,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))

_params_sk_silva_guillen_mose2 = SKSimpleParametersList(dict(zip(
    # from https://www.mdpi.com/2076-3417/6/10/284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.288,          0.089,          0.256,          "MoSe2",        0.720,
        -1.144,         -0.250,         -1.488,         -4.931,         -7.503,
        3.728,          -1.222,         -0.823,         0.215,          0.192,
        1.256,          -0.205,         1.256,          -0.205
    ]
)))

_params_sk_silva_guillen_ws2 = SKSimpleParametersList(dict(zip(
    # from https://www.mdpi.com/2076-3417/6/10/284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.153,          0.271,          0.057,          "WS2",          0.712,
        -1.155,         -0.650,         -2.279,         -3.864,         -7.327,
        7.911,          -1.220,         -1.328,         0.121,          0.422,
        1.178,          -0.273,         1.178,          -0.273
    ]
)))

_params_sk_silva_guillen_wse2 = SKSimpleParametersList(dict(zip(
    # from https://www.mdpi.com/2076-3417/6/10/284
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        3.260,          0.251,          0.439,          "WSe2",         0.722,
        -0.935,         -1.250,         -2.321,         -5.629,         -6.759,
        5.803,          -1.081,         -1.129,         0.094,          0.317,
        1.530,          -0.123,         1.530,          -0.123
    ]
)))

_params_sk_pearce = SKSimpleParametersList(dict(zip(
    # from https://doi.org/10.1103/PhysRevB.94.155416
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.319,          0.075,          0.052,          "MoS2",         0.687,
        2.12,           -0.46,          -1.41,          -5.38,          -3.69,
        -2.83,          0.67,           -0.24,          -0.62,          0.45,
        -0.42,          -1.32,          -0.42,          -1.32
    ]
)))

_params_sk_bieniek = SKSimpleParametersList(dict(zip(
    # from https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.085153
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3193,         0.075,          0.052,          "MoS2",         0.71,
        -0.03,          -0.03,          -0.03,          -3.36,          -4.78,
        -3.39,          1.1,            -1.1,           0.76,           0.27,
        1.19,           -0.83,          1.19,           -0.83
    ]
)))

_params_sk_abdi_mos2 = SKSimpleParametersList(dict(zip(
    # from https://iopscience.iop.org/article/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3166,         0.086,          0.052,          "MoS2",         0.710,
        -1.094,         -0.050,         -1.511,         -3.559,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))

_params_sk_abdi_wse2 = SKSimpleParametersList(dict(zip(
    # from https://iopscience.iop.org/article/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3260,         0.251,          0.439,          "WSe2",         0.722,
        -0.935,         -1.250,         -2.321,         -5.629,         -6.759,
        5.803,          -1.081,         -1.129,         0.094,          0.317,
        1.530,          -0.123,         1.530,          -0.205
    ]
)))

_params_sk_abdi_mose2 = SKSimpleParametersList(dict(zip(
    # from https://iopscience.iop.org/article/10.1149/2162-8777/abb28b
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.3298,         0.089,          0.256,          "MoSe2",        0.721,
        -1.144,         -0.250,         -1.488,         -4.931,         -7.503,
        3.728,          -1.222,         -0.823,         0.215,          0.192,
        1.256,          -0.205,         1.256,          -0.205,
    ]
)))

_params_sk_rostami = SKSimpleParametersList(dict(zip(
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_1_pds",      "v_1_pdp",      "v_2_dds",      "v_2_ddp",      "v_2_ddd",
        "v_0_pps",      "v_0_ppp",      "v_2_pps",      "v_2_ppp"
    ],
    [
        0.316,          0.075,          0.052,          "MoS2",         0.716,
        -1.094,         None,           -1.512,         -3.560,         -6.886,
        3.689,          -1.241,         -0.895,         0.252,          0.228,
        1.225,          -0.467,         1.225,          -0.467
    ]
)))
