def prepare_Hill1948_parameters(parameters):
    """Prepare the parameter list for the Hill1948 yield function."""
    # TODO: correct this; should be given in terms of "anisotropic yield stress ratios"
    return [
        parameters['F'],
        parameters['G'],
        parameters['H'],
        parameters['L'],
        parameters['M'],
        parameters['N'],
    ],


def prepare_Barlat_Yld91_parameters(parameters):
    """Prepare the parameter list for the Barlat_Yld91 yield function."""
    return [
        parameters['a'],
        parameters['b'],
        parameters['c'],
        parameters['f'],
        parameters['g'],
        parameters['h'],
    ]


def prepare_Barlat_Yld2004_18p_parameters(parameters):
    """Prepare the parameter list for the Barlat_Yld2004_18p yield function."""
    return [
        parameters['c_p_12'],
        parameters['c_p_13'],
        parameters['c_p_21'],
        parameters['c_p_23'],
        parameters['c_p_31'],
        parameters['c_p_32'],
        parameters['c_p_66'] * 0.5,
        parameters['c_p_55'] * 0.5,
        parameters['c_p_44'] * 0.5,
        parameters['c_dp_12'],
        parameters['c_dp_13'],
        parameters['c_dp_21'],
        parameters['c_dp_23'],
        parameters['c_dp_31'],
        parameters['c_dp_32'],
        parameters['c_dp_66'] * 0.5,
        parameters['c_dp_55'] * 0.5,
        parameters['c_dp_44'] * 0.5,
    ]


YIELD_FUNC_LOOKUP = {
    'Barlat_Yld91': prepare_Barlat_Yld91_parameters,
    'Barlat_Yld2004_18p': prepare_Barlat_Yld2004_18p_parameters,
}
