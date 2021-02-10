import numpy as np

def Generate_FE_sample(sample_size, Inhomogeneity_factor, L_groove, L_slope, 
                       Material_angle, Groove_angle,Elastic_modulus, Poisson_ratio,
                       Density, Barlat, Path_plastic_table):
    
    Plastic = np.loadtxt(Path_plastic_table, comments='%', delimiter=',')
       
    Sample_input_data = [sample_size, Inhomogeneity_factor, L_groove, L_slope, 
                       Material_angle, Groove_angle, Elastic_modulus, Poisson_ratio,
                       Density, Barlat, Plastic]
    
    return Sample_input_data
    
    
































