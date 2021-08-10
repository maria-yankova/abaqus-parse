from abaqus_parse import write_MK_mesh


def generate_MK_mesh(name_inp, FE_input_data):

    ######    Sample data   ######
    sample_size = FE_input_data['sample_size']
    # Groove geometry
    # thickness of the groove
    Inhomogeneity_factor = FE_input_data['inhomogeneity_factor']
    L_groove = FE_input_data['L_groove']
    L_slope = FE_input_data['L_slope']
    # Material angle
    Material_angle = FE_input_data['material_angle']  # (°)
    # Groove angle
    Groove_angle = 90 - FE_input_data['groove_angle']  # (°)
    # Material properties
    E = FE_input_data['elastic_modulus']
    mu = FE_input_data['poisson_ratio']
    rho = FE_input_data['density']

    law = FE_input_data['law']

    Plastic = FE_input_data['plastic']

    ######    FE data   ######
    mesh_size = FE_input_data['mesh_size']
    Element_type = FE_input_data['elem_type']
    Nb_el_thickness = FE_input_data['Nb_el_thickness']
    # Strain rate (/s)
    Eps_rate = FE_input_data['strain_rate']

    # Step time (put 0 => default Abaqus value)
    time_step = FE_input_data['total_time']  # Time period of the step (s)
    dt_i = FE_input_data['time_step'][0]  # Suggested initial time increment (s)
    dt_min = FE_input_data['time_step'][1]  # Minimum time increment allowed (s)
    dt_max = FE_input_data['time_step'][2]  # Maximum time increment allowed (s)
    # Bulk viscosity
    b1 = FE_input_data['bulk_parameters'][0]  # Linear bulk viscosity parameter
    b2 = FE_input_data['bulk_parameters'][1]  # Quadratic bulk viscosity parameter
    # Boundary conditions
    U_left = FE_input_data['displacment_BC'][0]  # (mm) along x
    U_right = FE_input_data['displacment_BC'][1]  # (mm) along x
    U_up = FE_input_data['displacment_BC'][2]  # (mm) along y
    U_bottom = FE_input_data['displacment_BC'][3]  # (mm) along y
    Max_plastic_strain = FE_input_data['max_plastic_strain']  # Eq. plastic strain PEEQ
    num_interval = FE_input_data.get('num_interval')

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                             Writing
    # #############################################################################
    write_MK_mesh.heading(name_inp, Inhomogeneity_factor,
                          Material_angle, Groove_angle, Eps_rate, law)

    write_MK_mesh.part(name_inp, sample_size, mesh_size, Inhomogeneity_factor,
                       L_groove, L_slope, Element_type, Material_angle, Groove_angle, Nb_el_thickness)

    write_MK_mesh.amplitude(name_inp, Eps_rate)

    write_MK_mesh.material(name_inp, E, mu, rho, Plastic, law)

    write_MK_mesh.step(name_inp, time_step, dt_i, dt_min, dt_max, b1, b2, U_left,
                       U_right, U_up, U_bottom, Max_plastic_strain, num_interval)
