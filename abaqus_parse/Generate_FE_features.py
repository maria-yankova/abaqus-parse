
def Generate_FE_features(mesh_size, bulk_parameters, elem_type, Strain_rate, total_time, Displacment_BC, Step_time):
   
    FE_input_data = [mesh_size, bulk_parameters, elem_type, Strain_rate, total_time, Displacment_BC, Step_time]
	
    return FE_input_data
    
