from abaqus_parse import write_MK_mesh

def Generate_MK_mesh(name_inp, Sample_input_data, FE_input_data):
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                             Inputs
    # #############################################################################

                       
    ######    Sample data   ###### 
    sample_size = Sample_input_data[0]    
    # Groove geometry
    Inhomogeneity_factor = Sample_input_data[1]  # thickness of the groove
    L_groove = Sample_input_data[2]
    L_slope = Sample_input_data[3]    
    # Material angle
    Material_angle = Sample_input_data[4] # (°)    
    # Groove angle
    Groove_angle = Sample_input_data[5] # (°) 
    # Material properties
    E = Sample_input_data[6]  
    mu = Sample_input_data[7]
    rho = Sample_input_data[8]
    power = Sample_input_data[9][0][0] 
    Barlat = Sample_input_data[9][1:]  
    Plastic = Sample_input_data[10]  

    
    ######    FE data   ###### 
    mesh_size = FE_input_data[0]
    Element_type = FE_input_data[2]   
    # Strain rate (/s)
    Eps_rate = FE_input_data[3]  
    # Step time (put 0 => default Abaqus value)
    time_step = FE_input_data[4]  # Time period of the step (s)
    dt_i = FE_input_data[6][0] # Suggested initial time increment (s)
    dt_min = FE_input_data[6][1]  # Minimum time increment allowed (s)
    dt_max = FE_input_data[6][2]  # Maximum time increment allowed (s)
    # Bulk viscosity
    b1 = FE_input_data[1][0]  # Linear bulk viscosity parameter
    b2 = FE_input_data[1][1]  # Quadratic bulk viscosity parameter  
    # Boundary conditions
    U_left = FE_input_data[5][0]  # (mm) along x
    U_right = FE_input_data[5][1]  # (mm) along x
    U_up = FE_input_data[5][2]  # (mm) along y
    U_bottom = FE_input_data[5][3] # (mm) along y
    
    
    
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #                             Writing
    # #############################################################################
    write_MK_mesh.heading(name_inp, Inhomogeneity_factor, Material_angle, Groove_angle, Eps_rate)
    
    write_MK_mesh.part(name_inp, sample_size, mesh_size, Inhomogeneity_factor,
                   L_groove, L_slope, Element_type, Material_angle, Groove_angle)
    
    write_MK_mesh.amplitude(name_inp, Eps_rate)
    
    write_MK_mesh.material(name_inp, E, mu, rho, Plastic, power, Barlat)
    
    write_MK_mesh.step(name_inp, time_step, dt_i, dt_min, dt_max, b1, b2, U_left,
                   U_right, U_up, U_bottom)
