import numpy as np


def generate_FE_input(sample_size, inhomogeneity_factor, L_groove, L_slope, material_angle, groove_angle, elastic_modulus, poisson_ratio, density, law, path_plastic_table, mesh_size, bulk_parameters, elem_type, strain_rate, total_time, displacment_BC, time_step, nb_el_thickness, max_plastic_strain):

    plastic = np.loadtxt(path_plastic_table, comments='%', delimiter=',')

    FE_input_data = [sample_size, inhomogeneity_factor, L_groove, L_slope,
                     material_angle, groove_angle, elastic_modulus, poisson_ratio,
                     density, law, plastic, mesh_size, bulk_parameters, elem_type,
                     strain_rate, total_time, displacment_BC, time_step,
                     nb_el_thickness, max_plastic_strain]

    return FE_input_data
