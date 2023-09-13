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
    # from https://link.aps.org/doi/10.1103/PhysRevB.92.195402
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

_params_sk_dias_mos2 = SKParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3166,         0.0806,         0.0536,         "MoS2",         0.710,
        -0.4939,        0.5624,         -0.2473,        -3.14115,       -4.0595,
        4.1989,         -1.89235,
        4.2398,         -1.2413,        2.2251,         -0.7614,
        -0.6717,        0.5706,         0.2729,         -0.0914,        -0.4619,
        0.0150,         0.0497,         0.8131,         -0.2763,
        0.0314,         0.0961,         -0.0305,        0.3723,         0.0014,
        0.0051,         0.0184,         -0.0395,        0.0092
    ]
)))

_params_sk_dias_mose2 = SKParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3288,         0.0806,         0.0820,         "MoSe2",        0.710,
        -0.1276,        0.3046,         -0.2724,        -3.8352,        -4.30195,
        4.30095,        -2.8093,
        3.4524,         -1.4295,        2.0197,         -0.6811,
        -0.6674,        0.5573,         0.0970,         1.2630,         -0.4857,
        0.01637,        0.0965,         0.9449,         -0.3039,
        0.0776,         0.0573,         -0.04778,       0.2372,         0.0249,
        0.0140,         0.0354,         -0.0293,        -0.0094
    ]
)))

_params_sk_dias_mote2 = SKParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3519,         0.0806,         0.1020,         "MoTe2",        0.710,
        -0.6630,        0.0491,         -0.2852,        -0.9084,        -1.8434,
        2.6799,         0.0678,
        2.2362,         -0.6279,        1.8294,         -0.5048,
        -0.4795,        -0.0934,        0.1656,         0.8198,         -0.2483,
        0.3267,         0.3033,         0.8459,         -0.4143,
        -0.1493,        -0.0627,        0.0360,         0.1169,         0.2683,
        -0.0617,        0.1002,         0.0114,         -0.0092
    ]
)))

_params_sk_dias_ws2 = SKParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.31532,        0.2754,         0.0536,         "WS2",          0.710,
        -0.3609,        0.8877,         -0.7364,        -3.52825,       -4.5926,
        4.415,          -1.97685,
        5.2769,         -1.2119,        2.4044,         -0.8115,
        -0.8942,        0.7347,         0.3417,         -0.3943,        -0.4069,
        -0.0142,        0.0036,         0.8415,         -0.2661,
        0.0508,         0.1278,         -0.0091,        0.1415,         0.0261,
        -0.0135,        -0.0191,        -0.0169,        0.0262
    ]
)))

_params_sk_dias_wse2 = SKParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.98.075202
    [
        "a",            "lamb_m",       "lamb_x",       "material",     "theta",
        "delta_0",      "delta_1",      "delta_2",      "delta_p",      "delta_z",
        "v_0_pps",      "v_0_ppp",
        "v_1_e_pds",    "v_1_e_pdp",    "v_1_o_pds",    "v_1_o_pdp",
        "v_2_e_dds",    "v_2_e_ddp",    "v_2_e_ddd",    "v_2_e_pps",    "v_2_e_ppp",
        "v_2_o_ddp",    "v_2_o_ddd",    "v_2_o_pps",    "v_2_o_ppp",
        "v_5_e_dds",    "v_5_e_ddp",    "v_5_e_ddd",    "v_5_e_pps",    "v_5_e_ppp",
        "v_5_o_ddp",    "v_5_o_ddd",    "v_5_o_pps",    "v_5_o_ppp"
    ],
    [
        0.3282,         0.2754,         0.0820,         "WSe2",         0.710,
        -0.5558,        0.6233,         -1.9340,        -2.20345,       -3.6177,
        3.1056,         -0.99385,
        5.1750,         -0.9139,        2.1733,         -0.7688,
        -0.8697,        0.6206,         0.3743,         0.1311,         -0.2475,
        -0.0469,        0.0923,         0.9703,         -0.2920,
        0.0443,         0.0912,         -0.0447,        0.1197,         0.1075,
        0.0096,         0.0140,         -0.0451,        0.0113
    ]
)))
