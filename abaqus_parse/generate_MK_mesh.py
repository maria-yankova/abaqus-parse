from abaqus_parse import write_MK_mesh


def generate_MK_mesh(name_inp, FE_input_data):

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                             Inputs
    # #############################################################################

    ######    Sample data   ######
    sample_size = FE_input_data[0]
    # Groove geometry
    Inhomogeneity_factor = FE_input_data[1]  # thickness of the groove
    L_groove = FE_input_data[2]
    L_slope = FE_input_data[3]
    # Material angle
    Material_angle = FE_input_data[4]  # (°)
    # Groove angle
    Groove_angle = 90 - FE_input_data[5]  # (°)
    # Material properties
    E = FE_input_data[6]
    mu = FE_input_data[7]
    rho = FE_input_data[8]
    law = FE_input_data[9][0][0]
    power = FE_input_data[9][1][0]
    Barlat = FE_input_data[9][2]
    Plastic = FE_input_data[10]

    ######    FE data   ######
    mesh_size = FE_input_data[11]
    Element_type = FE_input_data[13]
    Nb_el_thickness = FE_input_data[18]
    # Strain rate (/s)
    Eps_rate = FE_input_data[14]
    # Step time (put 0 => default Abaqus value)
    time_step = FE_input_data[15]  # Time period of the step (s)
    dt_i = FE_input_data[17][0]  # Suggested initial time increment (s)
    dt_min = FE_input_data[17][1]  # Minimum time increment allowed (s)
    dt_max = FE_input_data[17][2]  # Maximum time increment allowed (s)
    # Bulk viscosity
    b1 = FE_input_data[12][0]  # Linear bulk viscosity parameter
    b2 = FE_input_data[12][1]  # Quadratic bulk viscosity parameter
    # Boundary conditions
    U_left = FE_input_data[16][0]  # (mm) along x
    U_right = FE_input_data[16][1]  # (mm) along x
    U_up = FE_input_data[16][2]  # (mm) along y
    U_bottom = FE_input_data[16][3]  # (mm) along y
    Max_plastic_strain = FE_input_data[19]  # Eq. plastic strain PEEQ

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                             Writing
    # #############################################################################
    write_MK_mesh.heading(name_inp, Inhomogeneity_factor,
                          Material_angle, Groove_angle, Eps_rate, law)

    write_MK_mesh.part(name_inp, sample_size, mesh_size, Inhomogeneity_factor,
                       L_groove, L_slope, Element_type, Material_angle, Groove_angle, Nb_el_thickness)

    write_MK_mesh.amplitude(name_inp, Eps_rate)

    write_MK_mesh.material(name_inp, E, mu, rho, Plastic, power, Barlat, law)

    write_MK_mesh.step(name_inp, time_step, dt_i, dt_min, dt_max, b1, b2, U_left,
                       U_right, U_up, U_bottom, Max_plastic_strain)
