from .tmd_abstract_lattice import ParametersList


_params_sg_liu_2nn_mos2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.319,          0.073,          None,           "MoS2",
        1.046,          2.104,
        -0.184,         0.507,          -0.401,         0.057,          0.338,          0.218
    ]
)))

_params_sg_liu_2nn_mose2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3326,         0.091,          None,           "MoSe2",
        0.919,          2.065,
        -0.188,         0.456,          -0.317,         0.13,           0.29,           0.211
    ]
)))

_params_sg_liu_2nn_mote2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3557,         0.107,          None,           "MoTe2",
        0.605,          1.972,
        -0.169,         0.39,           -0.228,         0.252,          0.239,          0.207,
    ]
)))

_params_sg_liu_2nn_ws2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3191,         0.211,          None,           "WS2",
        1.13,           2.275,
        -0.206,         0.536,          -0.567,         -0.061,         0.384,          0.286
    ]
)))

_params_sg_liu_2nn_wse2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.3325,         0.228,          None,           "WSe2",
        0.943,          2.179,
        -0.207,         0.486,          -0.457,         0.034,          0.329,          0.263
    ]
)))

_params_sg_liu_2nn_wte2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
    ],
    [
        0.356,          0.237,          None,           "WTe2",
        0.606,          2.102,
        -0.175,         0.41,           -0.342,         0.19,           0.27,           0.233
    ]
)))

_params_sg_liu_6nn_mos2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3190,         0.073,          None,           "MoS2",
        0.683,          1.707,
        -0.146,         0.506,          0.114,          0.073,          0.162,          0.085,
        0.060,          0.273,          -0.034,         0.167,          -0.077,
        -0.038,         0.001,          -0.046,         -0.150,         -0.176,         0.266
    ]
)))

_params_sg_liu_6nn_mose2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3326,         0.091,          None,           "MoSe2",
        0.684,          1.546,
        -0.146,         0.432,          0.13,           0.075,          0.117,          0.144,
        0.039,          0.2413,         0.0174,         0.1559,         -0.0797,
        -0.042,         0.008,          -0.036,         -0.15,          -0.172,         0.272,
    ]
)))

_params_sg_liu_6nn_mote2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3557,         0.107,          None,           "MoTe2",
        0.588,          1.303,
        -0.226,         0.036,          0.234,          0.017,          0.098,          0.4,
        0.003,          0.0289,         0.0526,         0.1703,         0.1951,
        0.057,          0.187,          -0.103,         0.087,          -0.141,         -0.045,
    ]
)))

_params_sg_liu_6nn_ws2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3191,         0.211,          None,           "WS2",
        0.717,          1.916,
        -0.152,         0.59,           0.097,          0.016,          0.178,          0.047,
        0.069,          0.3014,         -0.0659,        0.1858,         -0.1236,
        -0.054,         0.002,          -0.045,         -0.163,         -0.206,         0.325,
    ]
)))

_params_sg_liu_6nn_wse2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.3325,         0.228,          None,           "WSe2",
        0.728,          1.655,
        -0.146,         0.507,          0.124,          0.015,          0.127,          0.117,
        0.036,          0.2702,         0.0007,         0.1739,         -0.1236,
        -0.061,         0.007,          -0.032,         -0.164,         -0.202,         0.329,
    ]
)))

_params_sg_liu_6nn_wte2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.88.085433
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e"
    ],
    [
        0.356,          0.237,          None,           "WTe2",
        0.697,          1.38,
        -0.109,         0.368,          0.164,          0.038,          0.093,          0.204,
        -0.015,         0.241,          0.110,          0.1306,         -0.1236,
        -0.066,         -0.013,         -0.011,         -0.132,         -0.177,         0.312,
    ]
)))

_params_sg_wu = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.91.075310
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_5_0_m_e",    "u_5_1_m_e",    "u_5_3_m_e",    "u_5_5_m_e",    "u_5_6_m_e",
        "u_5_0_m_o",    "u_5_2_m_o",
        "u_6_0_m_e",    "u_6_1_m_e",    "u_6_2_m_e",    "u_6_3_m_e",    "u_6_4_m_e",    "u_6_5_m_e",
        "u_6_0_m_o",    "u_6_1_m_o",    "u_6_2_m_o",
    ],
    [
        0.3190,         0.073,          None,           "MoS2",
        0.683,          1.707,          3.558,
        -0.146,         0.506,          0.114,          0.073,          0.162,          0.085,
        -0.189,         -0.024,         -0.117,
        0.060,          0.273,          -0.034,         0.167,          -0.077,
        -0.063,         0.025,
        -0.038,         0.001,          -0.046,          -0.150,         -0.176,         0.266,
        0.165,          0.140,          -0.122
    ]
)))

_params_sg_fang_mos2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.318,          0.0836,         0.0556,         "MoS2",
        -0.138,         0.0874,         1.0688,
        -1.9065,        -2.8949,        -1.2902,        -0.7755,
        1.4114,         -0.9402,        0.6517,         -0.8836,        -0.9535,
        -0.7883,        2.1584,         -1.379,
        -0.2979,        0.4096,         0.1145,         -0.5581,        -0.2487,        0.2747,
        -0.2069,        0.2562,         0.0323,
        0.9122,         0.0385,         0.1063,         0.0059,         0.0075,         -0.1916,
        0.8651,         0.0705,         -0.0995,        -0.1872,        -0.0679,        -0.1739,
        None,           -0.1498,        -0.2451,        -0.0686,        -0.2205
    ]
)))

_params_sg_fang_ws2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.318,          0.2874,         0.0556,         "WS2",
        -0.0393,        0.1984,         1.3754,
        -2.3461,        -3.3706,        -1.5534,        -1.1278,
        1.5629,         -0.9878,        0.6718,         -1.013,         -0.9491,
        -0.8855,        2.3121,         -1.4376,
        -0.3716,        0.4896,         0.1467,         -0.6892,        -0.303,         0.3537,
        -0.2011,        0.3106,         0.0263,
        0.9673,         0.1018,         0.1645,         0.0143,         -0.0315,        -0.2112,
        0.8726,         0.0989,         -0.1105,        -0.2187,        -0.0818,        -0.1749,
        None,           -0.1533,        -0.2736,        -0.0659,        -0.2618
    ]
)))

_params_sg_fang_mose2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.332,          0.0836,         0.247,          "MoSe2",
        -0.2297,        0.0149,         0.7819,
        -1.7806,        -2.9015,        -1.1726,        -0.6567,
        1.2677,         -0.8738,        0.5545,         -0.772,         -0.8578,
        -0.6946,        1.9415,         -1.3258,
        -0.2636,        0.352,          0.096,          -0.4734,        -0.2012,        0.2505,
        -0.146,         0.1912,         0.0177,
        0.9911,         0.0394,         0.1216,         -0.0036,        0.0047,         -0.2166,
        0.9638,         0.068,          -0.0755,        -0.1724,        -0.0735,        -0.2112,
        None,           -0.1553,        -0.2154,        -0.0691,        -0.2227
    ]
)))

_params_sg_fang_wse2 = ParametersList(dict(zip(
    # from https://link.aps.org/doi/10.1103/PhysRevB.92.205108
    [
        "a",            "lamb_m",       "lamb_x",       "material",
        "eps_0_m_e",    "eps_1_m_e",    "eps_0_m_o",
        "eps_0_x_e",    "eps_1_x_e",    "eps_0_x_o",    "eps_1_x_o",
        "u_1_0_m_e",    "u_1_1_m_e",    "u_1_2_m_e",    "u_1_3_m_e",    "u_1_4_m_e",
        "u_1_0_m_o",    "u_1_1_m_o",    "u_1_2_m_o",
        "u_2_0_m_e",    "u_2_1_m_e",    "u_2_2_m_e",    "u_2_3_m_e",    "u_2_4_m_e",    "u_2_5_m_e",
        "u_2_0_m_o",    "u_2_1_m_o",    "u_2_2_m_o",
        "u_2_0_x_e",    "u_2_1_x_e",    "u_2_2_x_e",    "u_2_3_x_e",    "u_2_4_x_e",    "u_2_5_x_e",
        "u_2_0_x_o",    "u_2_1_x_o",    "u_2_2_x_o",    "u_2_3_x_o",    "u_2_4_x_o",    "u_2_5_x_o",
        "u_3_0_m_e",    "u_3_1_m_e",    "u_3_2_m_e",    "u_3_3_m_e",    "u_3_4_m_e",

    ],
    [
        0.332,          0.2874,         0.247,          "WSe2",
        -0.1667,        0.0984,         1.0349,
        -2.182,         -3.3642,        -1.3937,        -0.9573,
        1.403,          -0.9044,        0.5711,         -0.8998,        -0.8548,
        -0.7744,        2.0858,         -1.4014,
        -0.333,         0.4233,         0.125,          -0.5837,        -0.2456,        0.319,
        -0.1395,        0.2321,         0.0129,
        1.047,          0.1027,         0.1857,         0.0029,         -0.0377,        -0.2399,
        0.9763,         0.092,          -0.0797,        -0.1985,        -0.0912,        -0.2171,
        None,           -0.1608,        -0.2424,        -0.0676,        -0.2618
    ]
)))
