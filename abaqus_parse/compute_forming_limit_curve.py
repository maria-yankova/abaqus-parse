import numpy as np

def compute_forming_limit_curve(all_model_responses):
    
    time = []
    S_mid_mises = []
    LE_mid_mises = []
    LE_mid_11 = []
    LE_mid_22 = []
    S_corner_mises = []
    LE_corner_mises = []
    LE_corner_11 = []
    LE_corner_22 = []
    
    dt = []
    first_derivative_mid = []
    first_derivative_corner = []
    
    threshold = 10
    
    minor_strain = []
    major_strain = []
    
    
    for i in range(len(all_model_responses)):
        time.append(all_model_responses[i][0])
        S_mid_mises.append(all_model_responses[i][1])
        LE_mid_mises.append(all_model_responses[i][2])
        LE_mid_11.append(all_model_responses[i][3])
        LE_mid_22.append(all_model_responses[i][4])
        S_corner_mises.append(all_model_responses[i][5])
        LE_corner_mises.append(all_model_responses[i][6])
        LE_corner_11.append(all_model_responses[i][7])
        LE_corner_22.append( all_model_responses[i][8])

        dt.append(np.gradient(time[i]))
        first_derivative_mid.append(np.gradient(LE_mid_mises[i])/dt[i])
        first_derivative_corner.append(np.gradient(LE_corner_mises[i])/dt[i])
        
        strain_rate_ratio = first_derivative_mid[i]/first_derivative_corner[i]       
					
        idx_necking = next(x for x, val in enumerate(strain_rate_ratio) if val > threshold)
        
        minor_strain.append(LE_corner_22[i][0:idx_necking])
        major_strain.append(LE_corner_11[i][0:idx_necking])
    
    forming_limit_curve = [minor_strain, major_strain]
    
    with open('FLC.dat', 'a') as inp_file:
        inp_file.write('** Data\n')
        for flc in forming_limit_curve:
            inp_file.write('** New_BC\n')
            for line in flc:
                str_line = ', '.join([' ' + str(round(x, 8)) for x in line])
                inp_file.write(str_line + '\n')

	
    return forming_limit_curve
    