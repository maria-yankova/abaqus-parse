import write_MK_mesh
import numpy as np

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                             Inputs
# #############################################################################

name_inp = 'Test.inp'

######    Sample data   ###### 
sample_size = [20,20,1]    
# Groove geometry
Inhomogeneity_factor = 0.9  # thickness of the groove
L_groove = 0.5
L_slope = 0.25 
# Material angle
Material_angle = 0 # (°)    
# Groove angle
Groove_angle = 90 - 30 # (°) 
# Material properties
E = 70000.0
mu = 0.33
rho = 2.7e-9
law = 'Barlat'
power = 5
Barlat =  [0.7087, 0.9074, 1.1376, 1.1526, 1.4851, 1.1637, 0.3039, 0.6731, 0.578, 0.9646, 1.0186, 1.4269, 1.1968, 0.5941, 0.6772, 0.6567,  0.2497, 0.41125] 
Plastic = np.loadtxt('Material_properties.txt', comments='%', delimiter=',')


######    FE data   ###### 
mesh_size = 0.3
Element_type = 'C3D8R'
# Strain rate (/s)
Eps_rate = 25
# Step time (put 0 => default Abaqus value)
time_step = 0.018 # Time period of the step (s)
dt_i = 0.0 # Suggested initial time increment (s)
dt_min = 0.0  # Minimum time increment allowed (s)
dt_max = 0.00001  # Maximum time increment allowed (s)
# Bulk viscosity
b1 = 0.1 # Linear bulk viscosity parameter
b2 = 0.1  # Quadratic bulk viscosity parameter  
# Boundary conditions
U_left = -0.5  # (mm) along x
U_right = 0.5 # (mm) along x
U_up = 0.5  # (mm) along y
U_bottom = -0.5 # (mm) along y




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                             Writing
# #############################################################################
write_MK_mesh.heading(name_inp, Inhomogeneity_factor, Material_angle, Groove_angle, Eps_rate, law)

write_MK_mesh.part(name_inp, sample_size, mesh_size, Inhomogeneity_factor,
               L_groove, L_slope, Element_type, Material_angle, Groove_angle)

write_MK_mesh.amplitude(name_inp, Eps_rate)

write_MK_mesh.material(name_inp, E, mu, rho, Plastic, power, Barlat, law)

write_MK_mesh.step(name_inp, time_step, dt_i, dt_min, dt_max, b1, b2, U_left,
               U_right, U_up, U_bottom)
