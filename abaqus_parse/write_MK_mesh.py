import numpy as np
import calfem.geometry as cfg
import calfem.mesh as cfm
# import calfem.vis as cfv
import warnings

from abaqus_parse.yield_functions import YIELD_FUNC_LOOKUP

def heading(file_name, Inhomogeneity_factor, Material_angle, Groove_angle, Eps_rate, law):
    
    with open(file_name, 'w') as inp_file:
        inp_file.write('*Heading\n')
        inp_file.write('**\n')
        inp_file.write('** Job name: Surfalex_18p\n')
        inp_file.write('**\n')
        inp_file.write('** FEATURES\n')
        inp_file.write('** Inhomogeneity factor: %s\n'%str(Inhomogeneity_factor))
        inp_file.write('** Material angle: %s (deg)\n'%str(Material_angle))
        inp_file.write('** Groove angle: %s (deg)\n'%str(Groove_angle))
        inp_file.write('** Strain rate: %s (/s)\n'%str(Eps_rate))
        inp_file.write('** Law: %s \n'%law)
        inp_file.write('**\n')
        inp_file.write(' units : mm,N,kg\n')
        inp_file.write('**\n')
        inp_file.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO \n')
        inp_file.write('**\n')


def part(file_name, sample_size, mesh_size, Inhomogeneity_factor, L_groove,
         L_slope, Element_type, Material_angle, Groove_angle, Nb_el_thickness):

    with open(file_name, 'a') as inp_file:
        inp_file.write('** PART\n')
        inp_file.write('*Part,name=MK_sample\n')
        inp_file.write('**\n')
        inp_file.write('*End part\n')
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('** ASSEMBLY\n*Assembly,name=Assembly\n')
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('** INSTANCE\n*Instance, name=MK_sample-1, part=MK_sample\n')
        inp_file.write('**\n')
        inp_file.write('**\n')


    ###################################
    ######         Nodes         ######
    ###################################
    Rotate = 0
        
    if Groove_angle>45:
        Groove_angle = 90 - Groove_angle
        Rotate = 1
        
    L_slope = L_slope/np.cos(Groove_angle*np.pi/180)
    L_groove = L_groove/np.cos(Groove_angle*np.pi/180)
    
    
    if Groove_angle == 45:
        warnings.warn("Difficulties to mesh the sample. Please consider change the groove angle.")
        
    else:    
        Groove_angle = Groove_angle*np.pi/180
        
        g = cfg.Geometry()
        g.point([0.0, sample_size[1]/2-np.tan(Groove_angle)*sample_size[0]/2 + L_slope+L_groove/2]) # point 0
        g.point([sample_size[0],  np.tan(Groove_angle)*sample_size[0]/2+sample_size[1]/2 + L_slope+L_groove/2 ]) # point 1
        g.point([sample_size[0], sample_size[1]]) # point 2
        g.point([0.0, sample_size[1]]) # point 3
    
        g.spline([0, 1]) # line 0
        g.spline([1, 2]) # line 1
        g.spline([2, 3]) # line 2
        g.spline([3, 0]) # line 3
        
        g.surface([0, 1, 2, 3])
            
        mesh = cfm.GmshMesh(g)
        
        mesh.elType = 3          # Degrees of freedom per node.
        mesh.dofsPerNode = 1     # Factor that changes element sizes.
        mesh.elSizeFactor = mesh_size   # Element size Factor
        
        Corner = 0  
        test_mesh = 0
        while Corner == 0:
            test_mesh += 1
            coords, edof, dofs, bdofs, elementmarkers = mesh.create()
            Nb_fois_el = np.zeros(len(dofs))
            # Trouver les bords
            for x in edof:
                Nb_fois_el[x-1] = Nb_fois_el[x-1]+1
            if np.sum(Nb_fois_el==1)==4:
                Corner = 1
            else:
                mesh.elSizeFactor = mesh.elSizeFactor*0.95
            if test_mesh>10:
                warnings.warn("Difficulties to mesh the sample. Please consider refining.")
                break
                
        
        Z_coord = np.ones(len(coords))*sample_size[2]
        
        # Trouve les coordonnées des noeuds du bord
        i1 = next(x for x, value in enumerate(coords[3:-1,0]) if value>sample_size[0]-1e-5)
        Ege_bottom_1 = coords[4:3+i1,:]
        Ege_bottom_1 = np.vstack([coords[0,:],Ege_bottom_1,coords[1,:]])
        # Trouve le numéro des noeuds du bords
        Num_up = np.arange(4,3+i1)
        Num_up = np.append([0],Num_up)
        Num_up = np.append(Num_up,[1])
    
        Num_bottom = Num_up + len(coords)
        
        # Noeuds supérieurs de la groove
        Coord_groove_up = Ege_bottom_1.copy()
        Coord_groove_up[:,1] = Coord_groove_up[:,1]-L_slope
        
        Z_groove_up = np.ones(len(Coord_groove_up))*((sample_size[2]-Inhomogeneity_factor)/2+Inhomogeneity_factor)
        
        # Noeuds inférieurs de la groove
        Coord_groove_bottom = Coord_groove_up.copy()
        Coord_groove_bottom[:,1] = Coord_groove_bottom[:,1]-L_groove
    
        Z_groove_bottom = np.ones(len(Coord_groove_bottom))*((sample_size[2]-Inhomogeneity_factor)/2+Inhomogeneity_factor)

        # Deuxième partie de l'éprouvette
        coords_2 = coords.copy()
        coords_2[:,0] = sample_size[0]-coords[:,0]
        coords_2[:,1] = sample_size[1]-coords[:,1]
       
        Z_coord_2 = Z_coord.copy()
        
        # On rassemble tous les noeuds d'une face
        coord_all = coords.copy()
        coord_all = np.vstack([coords,coords_2,Coord_groove_up,Coord_groove_bottom])
        
        Z_coord_all = np.concatenate([Z_coord,Z_coord_2,Z_groove_up,Z_groove_bottom])
    
        edof_2 = edof.copy()
        edof_2 = edof_2+np.max(edof_2)
        edof_all = edof.copy()
        edof_all = np.vstack([edof,edof_2])
        
        connectivite = np.zeros((len(Num_up)-1,4))
        
        # Rétrecissement 1
        for i in range(len(Num_up)-1):
            connectivite[i,:] = [len(coords)+len(coords_2)+i,len(coords)+len(coords_2)+i+1,Num_up[i+1],Num_up[i]]
            
        # Groove    
        connectivite2 = np.zeros((len(Num_up)-1,4))
        for i in range(len(Num_up)-1):
            connectivite2[i,:] = [len(coords)+len(coords_2)+len(Coord_groove_up)+i,
                                  len(coords)+len(coords_2)+len(Coord_groove_up)+i+1,len(coords)+len(coords_2)+i+1,len(coords)+len(coords_2)+i]    
        
        # Rétrecissement 2
        connectivite3 = np.zeros((len(Num_up)-1,4))
        for i in range(len(Num_up)-1):
            connectivite3[i,:] = [Num_bottom[-i-1],Num_bottom[-i-2],len(coords)+len(coords_2)+len(Coord_groove_up)+i+1,len(coords)+len(coords_2)+len(Coord_groove_up)+i]
            
       
        edof_all = np.vstack([edof_all,connectivite+1,connectivite2+1,connectivite3+1])
                
        if Rotate==1:
            C = coord_all.copy()
            C[:,0] = -coord_all[:,1] + sample_size[1]
            C[:,1] = coord_all[:,0]
            C[:,0] = sample_size[1]-C[:,0]   
            coord_all[:,0] = C[:,0] 
            coord_all[:,1] = C[:,1]

    
    # cfv.figure()
    # cfv.drawMesh(
    #     coords=coord_all,
    #     edof=edof_all,
    #     dofs_per_node=mesh.dofsPerNode,
    #     el_type=mesh.elType,
    #     filled=True,
    #     title="Mesh"   
    #         )   
    
    coord_all_th = []
    Z_coord_th = Z_coord_all.copy()
    X_coord_th = coord_all[:,0]
    Y_coord_th = coord_all[:,1]
    
    Nb_el_thickness = int(Nb_el_thickness)
    Coordinate = np.zeros(((Nb_el_thickness+1)*len(coord_all), 4))   
    Coordinate[:,0] = np.array(range((Nb_el_thickness+1)*len(coord_all)))+1
    

    
    for i in range(Nb_el_thickness):
        Z_coord_all_surface = Z_coord_all.copy()
        Z_coord_all_surface = abs(Z_coord_all-sample_size[2])
        
        coord_all_th.append(coord_all.copy())
        
        Z_coord_th  = np.concatenate([Z_coord_th, Z_coord_all * (Nb_el_thickness - (i+1))/Nb_el_thickness + Z_coord_all_surface * (i+1)/Nb_el_thickness])
        X_coord_th = np.concatenate([X_coord_th,coord_all[:,0]])
        Y_coord_th = np.concatenate([Y_coord_th,coord_all[:,1]])

                          
    Coordinate[:,1] = X_coord_th
    Coordinate[:,2] = Y_coord_th
    Coordinate[:,3] = Z_coord_th
    
    
        

    if Rotate==1:
        Coordinate[:,3] = sample_size[2]-Coordinate[:,3]
            

    with open(file_name, 'a') as inp_file:
        inp_file.write('*Node, nset=Sheet_nodes\n')
        for line in Coordinate:
            inp_file.write('%i, %.6f, %.6f, %.6f\n'%(line[0],line[1],line[2],line[3]))
            
            

    ###################################
    ######       Elements        ######
    ###################################  
    Num_el = []
    el_node =  np.zeros((8))

    Num = 0
    for el in range(Nb_el_thickness):
        for i in range(len(edof_all)):
            Num += 1
            el_node[0] = edof_all[i,3] + len(coord_all)*el
            el_node[1] = edof_all[i,2] + len(coord_all)*el
            el_node[2] = edof_all[i,1] + len(coord_all)*el
            el_node[3] = edof_all[i,0] + len(coord_all)*el
            el_node[4] = edof_all[i,3] + len(coord_all)*(el+1)
            el_node[5] = edof_all[i,2] + len(coord_all)*(el+1)
            el_node[6] = edof_all[i,1] + len(coord_all)*(el+1)
            el_node[7] = edof_all[i,0] + len(coord_all)*(el+1)
            Num_el.append([Num] + list(el_node))
  
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('*Element, type=%s, elset=Sheet_Elements\n'%(Element_type))
        for line in Num_el:
            inp_file.write('%i, %i, %i, %i, %i, %i, %i, %i, %i\n'%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]))


    ##########################################
    ## Orientation - Section - End instance ##
    ##########################################
    Mat_Angle = Material_angle*np.pi/180
    Rot = [np.cos(Mat_Angle), -np.sin(Mat_Angle), 0.0, np.sin(Mat_Angle), np.cos(Mat_Angle), 0.0]
    str_line = ', '.join([' ' + str(round(x,5)) for x in Rot])

    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('*Orientation, name=Ori\n')
        inp_file.write(str_line + '\n')
        inp_file.write('3, 0.\n')
        inp_file.write('**\n')
        inp_file.write('** SECTION\n') 
        inp_file.write('*Solid Section, elset=Sheet_Elements, orientation=Ori, material=Material-1\n')
        inp_file.write(',\n')
        inp_file.write('*End Instance\n')


    ###################################
    ####           Sets           #####
    ###################################


    Left_nodes  = np.where(Coordinate[:,1]<1e-4)[0]+1
    Bottom_nodes = np.where(Coordinate[:,2]<1e-4)[0]+1
    Upper_nodes  = np.where(Coordinate[:,2]>(sample_size[1]-1e-4))[0]+1
    Right_nodes   = np.where(Coordinate[:,1]>(sample_size)[0]-1e-4)[0]+1
    
    Mid = abs(Coordinate[:,1]-sample_size[0]/2)+abs(Coordinate[:,2]-sample_size[1]/2)
    Mid = np.argmin(Mid) + 1 
    i = np.where(connectivite2==Mid)[0][0]
    idx_el_mid = len(edof_all) - len(connectivite3)-len(connectivite2)+i
    
    Corner =  1 
    idx_el_corner = Corner
    
    
    ### Bottom nodes
    Node_set = Bottom_nodes  
    nodetowrite = []
    [div,rest] = np.divmod(len(Node_set),16)   
    for i in range(div):
        nodetowrite.append(Node_set[i*16:i*16+16])      
    if rest!=0:
        nodetowrite.append(Node_set[div*16:div*16+rest])
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Node_Set : Bottom_nodes \n')
        inp_file.write('*nset, nset=Bottom_nodes, instance=MK_sample-1\n')
        for line in nodetowrite:
            str_line = ', '.join([' ' + str(x) for x in line])
            inp_file.write(str_line + '\n')


    ### Upper nodes
    Node_set = Upper_nodes  
    nodetowrite = []
    [div,rest] = np.divmod(len(Node_set),16)   
    for i in range(div):
        nodetowrite.append(Node_set[i*16:i*16+16])     
    if rest!=0:
        nodetowrite.append(Node_set[div*16:div*16+rest])
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Node_Set : Upper_nodes \n')
        inp_file.write('*nset, nset=Upper_nodes, instance=MK_sample-1\n')
        for line in nodetowrite:
            str_line = ', '.join([' ' + str(x) for x in line])
            inp_file.write(str_line + '\n')


    ### Left nodes
    Node_set = Left_nodes  
    nodetowrite = []
    [div,rest] = np.divmod(len(Node_set),16)   
    for i in range(div):
        nodetowrite.append(Node_set[i*16:i*16+16])  
    if rest!=0:
        nodetowrite.append(Node_set[div*16:div*16+rest])
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Node_Set : Left_nodes \n')
        inp_file.write('*nset, nset=Left_nodes, instance=MK_sample-1\n')
        for line in nodetowrite:
            str_line = ', '.join([' ' + str(x) for x in line])
            inp_file.write(str_line + '\n')


    ### Right nodes
    Node_set = Right_nodes 
    nodetowrite = []
    [div,rest] = np.divmod(len(Node_set),16)
    for i in range(div):
        nodetowrite.append(Node_set[i*16:i*16+16])
    if rest!=0:
        nodetowrite.append(Node_set[div*16:div*16+rest])
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Node_Set : Right_nodes \n')
        inp_file.write('*nset, nset=Right_nodes, instance=MK_sample-1\n')
        for line in nodetowrite:
            str_line = ', '.join([' ' + str(x) for x in line])
            inp_file.write(str_line + '\n')
    
    ### Middle element
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Element_Set : Middle_element \n')
        inp_file.write('*elset, elset=Middle_element, instance=MK_sample-1\n')
        inp_file.write('%i\n'%(idx_el_mid))
        
    ### Corner element
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('**Element_Set : Corner_element \n')
        inp_file.write('*elset, elset=Corner_element, instance=MK_sample-1\n')
        inp_file.write('%i\n'%(idx_el_corner))

    ###################################
    ####       End Assembly       #####
    ###################################
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('*End Assembly\n')
        
    

def amplitude(file_name, Eps_rate):

    Strain_Time_Table = [0.0, 0.0, 1.0/Eps_rate*0.01, 0.201003342, 1.0/Eps_rate*0.02,
        0.404026801, 1.0/Eps_rate*0.03, 0.609090679, 1.0/Eps_rate*0.04, 0.816215484, 1.0/Eps_rate*0.05,
        1.025421928, 1.0/Eps_rate*0.06, 1.236730931, 1.0/Eps_rate*0.07, 1.450163625, 1.0/Eps_rate*0.08,
        1.665741353, 1.0/Eps_rate*0.09, 1.883485674, 1.0/Eps_rate*0.1, 2.103418362, 1.0/Eps_rate*0.11,
        2.325561409, 1.0/Eps_rate*0.12, 2.549937032, 1.0/Eps_rate*0.13, 2.776567666, 1.0/Eps_rate*0.14,
        3.005475977, 1.0/Eps_rate*0.15, 3.236684855, 1.0/Eps_rate*0.16, 3.47021742, 1.0/Eps_rate*0.17,
        3.706097026, 1.0/Eps_rate*0.18, 3.944347262, 1.0/Eps_rate*0.19, 4.184991953, 1.0/Eps_rate*0.2,
        4.428055163, 1.0/Eps_rate*0.21, 4.673561199, 1.0/Eps_rate*0.22, 4.921534612, 1.0/Eps_rate*0.23,
        5.172000199, 1.0/Eps_rate*0.24, 5.424983006, 1.0/Eps_rate*0.25, 5.680508334, 1.0/Eps_rate*0.26,
        5.938601733, 1.0/Eps_rate*0.27, 6.199289015, 1.0/Eps_rate*0.28, 6.462596247, 1.0/Eps_rate*0.29,
        6.728549761, 1.0/Eps_rate*0.3, 6.997176152, 1.0/Eps_rate*0.31, 7.268502283, 1.0/Eps_rate*0.32,
        7.542555287, 1.0/Eps_rate*0.33, 7.819362569, 1.0/Eps_rate*0.34, 8.098951811, 1.0/Eps_rate*0.35,
        8.381350972, 1.0/Eps_rate*0.36, 8.666588291, 1.0/Eps_rate*0.37, 8.954692293, 1.0/Eps_rate*0.38,
        9.245691789, 1.0/Eps_rate*0.39, 9.539615878, 1.0/Eps_rate*0.4, 9.836493953, 1.0/Eps_rate*0.41,
        10.1363557, 1.0/Eps_rate*0.42, 10.43923111, 1.0/Eps_rate*0.43, 10.74515047, 1.0/Eps_rate*0.44,
        11.05414437, 1.0/Eps_rate*0.45, 11.36624371, 1.0/Eps_rate*0.46, 11.6814797, 1.0/Eps_rate*0.47,
        11.99988386, 1.0/Eps_rate*0.48, 12.32148804, 1.0/Eps_rate*0.49, 12.6463244, 1.0/Eps_rate*0.5,
        12.97442541, 1.0/Eps_rate*0.51, 13.3058239, 1.0/Eps_rate*0.52, 13.64055299, 1.0/Eps_rate*0.53,
        13.97864617, 1.0/Eps_rate*0.54, 14.32013724, 1.0/Eps_rate*0.55, 14.66506036, 1.0/Eps_rate*0.56,
        15.01345001, 1.0/Eps_rate*0.57, 15.36534103, 1.0/Eps_rate*0.58, 15.72076862, 1.0/Eps_rate*0.59,
        16.07976831, 1.0/Eps_rate*0.6, 16.44237601, 1.0/Eps_rate*0.61, 16.80862798, 1.0/Eps_rate*0.62,
        17.17856084, 1.0/Eps_rate*0.63, 17.55221159, 1.0/Eps_rate*0.64, 17.92961759, 1.0/Eps_rate*0.65,
        18.31081658, 1.0/Eps_rate*0.66, 18.69584669, 1.0/Eps_rate*0.67, 19.08474641, 1.0/Eps_rate*0.68,
        19.47755464, 1.0/Eps_rate*0.69, 19.87431066, 1.0/Eps_rate*0.7, 20.27505415, 1.0/Eps_rate*0.71,
        20.67982517, 1.0/Eps_rate*0.72, 21.08866421, 1.0/Eps_rate*0.73, 21.50161215, 1.0/Eps_rate*0.74,
        21.91871029, 1.0/Eps_rate*0.75, 22.34000033, 1.0/Eps_rate*0.76, 22.76552441, 1.0/Eps_rate*0.77,
        23.19532508, 1.0/Eps_rate*0.78, 23.62944531, 1.0/Eps_rate*0.79, 24.06792853, 1.0/Eps_rate*0.8,
        24.51081857, 1.0/Eps_rate*0.81, 24.95815973, 1.0/Eps_rate*0.82, 25.40999675, 1.0/Eps_rate*0.83,
        25.86637481, 1.0/Eps_rate*0.84, 26.32733954, 1.0/Eps_rate*0.85, 26.79293704, 1.0/Eps_rate*0.86,
        27.26321387, 1.0/Eps_rate*0.87, 27.73821707, 1.0/Eps_rate*0.88, 28.21799413, 1.0/Eps_rate*0.89,
        28.70259303, 1.0/Eps_rate*0.9, 29.19206222, 1.0/Eps_rate*0.91, 29.68645067, 1.0/Eps_rate*0.92,
        30.1858078, 1.0/Eps_rate*0.93, 30.69018355, 1.0/Eps_rate*0.94, 31.19962837, 1.0/Eps_rate*0.95,
        31.71419319, 1.0/Eps_rate*0.96, 32.23392947, 1.0/Eps_rate*0.97, 32.75888919, 1.0/Eps_rate*0.98,
        33.28912484, 1.0/Eps_rate*0.99, 33.82468945, 1.0/Eps_rate*1.0, 34.36563657, 1.0/Eps_rate*1.01,
        34.9120203, 1.0/Eps_rate*1.02, 35.46389528, 1.0/Eps_rate*1.03, 36.02131669, 1.0/Eps_rate*1.04,
        36.58434029, 1.0/Eps_rate*1.05, 37.15302236, 1.0/Eps_rate*1.06, 37.72741979, 1.0/Eps_rate*1.07,
        38.30759, 1.0/Eps_rate*1.08, 38.89359102, 1.0/Eps_rate*1.09, 39.48548145, 1.0/Eps_rate*1.1,
        40.08332048, 1.0/Eps_rate*1.11, 40.68716789, 1.0/Eps_rate*1.12, 41.29708407, 1.0/Eps_rate*1.13,
        41.91313, 1.0/Eps_rate*1.14, 42.5353673, 1.0/Eps_rate*1.15, 43.16385819, 1.0/Eps_rate*1.16,
        43.79866552, 1.0/Eps_rate*1.17, 44.43985277, 1.0/Eps_rate*1.18, 45.08748406, 1.0/Eps_rate*1.19,
        45.74162415, 1.0/Eps_rate*1.2, 46.40233845, 1.0/Eps_rate*1.21, 47.06969305, 1.0/Eps_rate*1.22,
        47.74375467, 1.0/Eps_rate*1.23, 48.42459073, 1.0/Eps_rate*1.24, 49.1122693, 1.0/Eps_rate*1.25,
        49.80685915, 1.0/Eps_rate*1.26, 50.50842975, 1.0/Eps_rate*1.27, 51.21705125, 1.0/Eps_rate*1.28,
        51.93279451, 1.0/Eps_rate*1.29, 52.65573112, 1.0/Eps_rate*1.3, 53.38593335, 1.0/Eps_rate*1.31,
        54.12347424, 1.0/Eps_rate*1.32, 54.86842755, 1.0/Eps_rate*1.33, 55.62086775, 1.0/Eps_rate*1.34,
        56.38087011, 1.0/Eps_rate*1.35, 57.14851061, 1.0/Eps_rate*1.36, 57.92386604, 1.0/Eps_rate*1.37,
        58.70701391, 1.0/Eps_rate*1.38, 59.49803255, 1.0/Eps_rate*1.39, 60.29700106, 1.0/Eps_rate*1.4,
        61.10399934, 1.0/Eps_rate*1.41, 61.91910808, 1.0/Eps_rate*1.42, 62.74240881, 1.0/Eps_rate*1.43,
        63.57398384, 1.0/Eps_rate*1.44, 64.41391634, 1.0/Eps_rate*1.45, 65.2622903, 1.0/Eps_rate*1.46,
        66.11919057, 1.0/Eps_rate*1.47, 66.98470282, 1.0/Eps_rate*1.48, 67.85891362, 1.0/Eps_rate*1.49,
        68.74191038, 1.0/Eps_rate*1.5, 69.63378141, 1.0/Eps_rate*1.51, 70.53461589, 1.0/Eps_rate*1.52,
        71.4445039, 1.0/Eps_rate*1.53, 72.36353645, 1.0/Eps_rate*1.54, 73.29180542, 1.0/Eps_rate*1.55,
        74.22940365, 1.0/Eps_rate*1.56, 75.1764249, 1.0/Eps_rate*1.57, 76.13296388, 1.0/Eps_rate*1.58,
        77.09911622, 1.0/Eps_rate*1.59, 78.07497857, 1.0/Eps_rate*1.6, 79.06064849, 1.0/Eps_rate*1.61,
        80.05622456, 1.0/Eps_rate*1.62, 81.06180633, 1.0/Eps_rate*1.63, 82.07749437, 1.0/Eps_rate*1.64,
        83.10339024, 1.0/Eps_rate*1.65, 84.13959654, 1.0/Eps_rate*1.66, 85.18621689, 1.0/Eps_rate*1.67,
        86.24335594, 1.0/Eps_rate*1.68, 87.31111942, 1.0/Eps_rate*1.69, 88.3896141, 1.0/Eps_rate*1.7,
        89.47894783, 1.0/Eps_rate*1.71, 90.57922955, 1.0/Eps_rate*1.72, 91.69056929, 1.0/Eps_rate*1.73,
        92.81307817, 1.0/Eps_rate*1.74, 93.94686845, 1.0/Eps_rate*1.75, 95.09205352, 1.0/Eps_rate*1.76,
        96.24874789, 1.0/Eps_rate*1.77, 97.41706723, 1.0/Eps_rate*1.78, 98.59712837, 1.0/Eps_rate*1.79,
        99.78904933, 1.0/Eps_rate*1.8, 100.9929493, 1.0/Eps_rate*1.81, 102.2089486, 1.0/Eps_rate*1.82,
        103.437169, 1.0/Eps_rate*1.83, 104.6777332, 1.0/Eps_rate*1.84, 105.9307652, 1.0/Eps_rate*1.85,
        107.1963905, 1.0/Eps_rate*1.86, 108.4747354, 1.0/Eps_rate*1.87, 109.765928, 1.0/Eps_rate*1.88,
        111.0700972, 1.0/Eps_rate*1.89, 112.3873736, 1.0/Eps_rate*1.9, 113.7178888, 1.0/Eps_rate*1.91,
        115.061776, 1.0/Eps_rate*1.92, 116.4191694, 1.0/Eps_rate*1.93, 117.7902048, 1.0/Eps_rate*1.94,
        119.1750194, 1.0/Eps_rate*1.95, 120.5737516, 1.0/Eps_rate*1.96, 121.9865413, 1.0/Eps_rate*1.97,
        123.4135298, 1.0/Eps_rate*1.98, 124.8548597, 1.0/Eps_rate*1.99, 126.3106752, 1.0/Eps_rate*2.0,
        127.781122, 1.0/Eps_rate*2.01, 129.2663469, 1.0/Eps_rate*2.02, 130.7664987, 1.0/Eps_rate*2.03,
        132.2817272, 1.0/Eps_rate*2.04, 133.812184, 1.0/Eps_rate*2.05, 135.3580221, 1.0/Eps_rate*2.06,
        136.9193962, 1.0/Eps_rate*2.07, 138.4964624, 1.0/Eps_rate*2.08, 140.0893783, 1.0/Eps_rate*2.09,
        141.6983033, 1.0/Eps_rate*2.1, 143.3233983, 1.0/Eps_rate*2.11, 144.9648257, 1.0/Eps_rate*2.12,
        146.6227498, 1.0/Eps_rate*2.13, 148.2973362, 1.0/Eps_rate*2.14, 149.9887526, 1.0/Eps_rate*2.15,
        151.6971679, 1.0/Eps_rate*2.16, 153.4227532, 1.0/Eps_rate*2.17, 155.1656808, 1.0/Eps_rate*2.18,
        156.9261252, 1.0/Eps_rate*2.19, 158.7042623, 1.0/Eps_rate*2.2, 160.50027, 1.0/Eps_rate*2.21,
        162.3143279, 1.0/Eps_rate*2.22, 164.1466173, 1.0/Eps_rate*2.23, 165.9973216, 1.0/Eps_rate*2.24,
        167.8666257, 1.0/Eps_rate*2.25, 169.7547167, 1.0/Eps_rate*2.26, 171.6617833, 1.0/Eps_rate*2.27,
        173.5880163, 1.0/Eps_rate*2.28, 175.5336082, 1.0/Eps_rate*2.29, 177.4987536, 1.0/Eps_rate*2.3,
        179.4836491, 1.0/Eps_rate*2.31, 181.4884931, 1.0/Eps_rate*2.32, 183.5134861, 1.0/Eps_rate*2.33,
        185.5588307, 1.0/Eps_rate*2.34, 187.6247313, 1.0/Eps_rate*2.35, 189.7113945, 1.0/Eps_rate*2.36,
        191.819029, 1.0/Eps_rate*2.37, 193.9478457, 1.0/Eps_rate*2.38, 196.0980573, 1.0/Eps_rate*2.39,
        198.2698789, 1.0/Eps_rate*2.4, 200.4635276, 1.0/Eps_rate*2.41, 202.6792229, 1.0/Eps_rate*2.42,
        204.9171863, 1.0/Eps_rate*2.43, 207.1776416, 1.0/Eps_rate*2.44, 209.4608149, 1.0/Eps_rate*2.45,
        211.7669344, 1.0/Eps_rate*2.46, 214.0962308, 1.0/Eps_rate*2.47, 216.448937, 1.0/Eps_rate*2.48,
        218.8252884, 1.0/Eps_rate*2.49, 221.2255224, 1.0/Eps_rate*2.5, 223.6498792, 1.0/Eps_rate*2.51,
        226.0986012, 1.0/Eps_rate*2.52, 228.5719333, 1.0/Eps_rate*2.53, 231.0701227, 1.0/Eps_rate*2.54,
        233.5934194, 1.0/Eps_rate*2.55, 236.1420757, 1.0/Eps_rate*2.56, 238.7163463, 1.0/Eps_rate*2.57,
        241.3164888, 1.0/Eps_rate*2.58, 243.9427632, 1.0/Eps_rate*2.59, 246.5954321, 1.0/Eps_rate*2.6,
        249.2747607, 1.0/Eps_rate*2.61, 251.981017, 1.0/Eps_rate*2.62, 254.7144717, 1.0/Eps_rate*2.63,
        257.475398, 1.0/Eps_rate*2.64, 260.2640722, 1.0/Eps_rate*2.65, 263.0807729, 1.0/Eps_rate*2.66,
        265.925782, 1.0/Eps_rate*2.67, 268.7993839, 1.0/Eps_rate*2.68, 271.7018659, 1.0/Eps_rate*2.69,
        274.6335184, 1.0/Eps_rate*2.7, 277.5946345, 1.0/Eps_rate*2.71, 280.5855103, 1.0/Eps_rate*2.72,
        283.6064449, 1.0/Eps_rate*2.73, 286.6577404, 1.0/Eps_rate*2.74, 289.7397019, 1.0/Eps_rate*2.75,
        292.8526377, 1.0/Eps_rate*2.76, 295.996859, 1.0/Eps_rate*2.77, 299.1726802, 1.0/Eps_rate*2.78,
        302.380419, 1.0/Eps_rate*2.79, 305.620396, 1.0/Eps_rate*2.8, 308.8929354, 1.0/Eps_rate*2.81,
        312.1983644, 1.0/Eps_rate*2.82, 315.5370134, 1.0/Eps_rate*2.83, 318.9092165, 1.0/Eps_rate*2.84,
        322.3153107, 1.0/Eps_rate*2.85, 325.7556368, 1.0/Eps_rate*2.86, 329.2305387, 1.0/Eps_rate*2.87,
        332.740364, 1.0/Eps_rate*2.88, 336.2854636, 1.0/Eps_rate*2.89, 339.866192, 1.0/Eps_rate*2.9,
        343.4829074, 1.0/Eps_rate*2.91, 347.1359713, 1.0/Eps_rate*2.92, 350.8257492, 1.0/Eps_rate*2.93,
        354.5526099, 1.0/Eps_rate*2.94, 358.3169262, 1.0/Eps_rate*2.95, 362.1190746, 1.0/Eps_rate*2.96,
        365.9594351, 1.0/Eps_rate*2.97, 369.8383919, 1.0/Eps_rate*2.98, 373.7563329, 1.0/Eps_rate*2.99,
        377.7136498, 1.0/Eps_rate*3.0, 381.7107385, 1.0/Eps_rate*3.01, 385.7479985, 1.0/Eps_rate*3.02,
        389.8258337, 1.0/Eps_rate*3.03, 393.9446518, 1.0/Eps_rate*3.04, 398.1048647, 1.0/Eps_rate*3.05,
        402.3068885, 1.0/Eps_rate*3.06, 406.5511432, 1.0/Eps_rate*3.07, 410.8380535, 1.0/Eps_rate*3.08,
        415.1680479, 1.0/Eps_rate*3.09, 419.5415595, 1.0/Eps_rate*3.1, 423.9590256, 1.0/Eps_rate*3.11,
        428.420888, 1.0/Eps_rate*3.12, 432.9275929, 1.0/Eps_rate*3.13, 437.4795908, 1.0/Eps_rate*3.14,
        442.0773372, 1.0/Eps_rate*3.15, 446.7212916, 1.0/Eps_rate*3.16, 451.4119186, 1.0/Eps_rate*3.17,
        456.1496871, 1.0/Eps_rate*3.18, 460.935071, 1.0/Eps_rate*3.19, 465.7685489, 1.0/Eps_rate*3.2,
        470.6506039, 1.0/Eps_rate*3.21, 475.5817245, 1.0/Eps_rate*3.22, 480.5624036, 1.0/Eps_rate*3.23,
        485.5931394, 1.0/Eps_rate*3.24, 490.6744349, 1.0/Eps_rate*3.25, 495.8067983, 1.0/Eps_rate*3.26,
        500.9907429, 1.0/Eps_rate*3.27, 506.2267869, 1.0/Eps_rate*3.28, 511.515454, 1.0/Eps_rate*3.29,
        516.8572731, 1.0/Eps_rate*3.3, 522.2527784, 1.0/Eps_rate*3.31, 527.7025094, 1.0/Eps_rate*3.32,
        533.2070112, 1.0/Eps_rate*3.33, 538.7668341, 1.0/Eps_rate*3.34, 544.3825341, 1.0/Eps_rate*3.35,
        550.0546729, 1.0/Eps_rate*3.36, 555.7838176, 1.0/Eps_rate*3.37, 561.5705412, 1.0/Eps_rate*3.38,
        567.4154223, 1.0/Eps_rate*3.39, 573.3190454, 1.0/Eps_rate*3.4, 579.2820009, 1.0/Eps_rate*3.41,
        585.3048852, 1.0/Eps_rate*3.42, 591.3883004, 1.0/Eps_rate*3.43, 597.532855, 1.0/Eps_rate*3.44,
        603.7391634, 1.0/Eps_rate*3.45, 610.0078462, 1.0/Eps_rate*3.46, 616.3395303, 1.0/Eps_rate*3.47,
        622.7348489, 1.0/Eps_rate*3.48, 629.1944415, 1.0/Eps_rate*3.49, 635.7189541, 1.0/Eps_rate*3.5,
        642.3090392, 1.0/Eps_rate*3.51, 648.9653557, 1.0/Eps_rate*3.52, 655.6885693, 1.0/Eps_rate*3.53,
        662.4793523, 1.0/Eps_rate*3.54, 669.3383838, 1.0/Eps_rate*3.55, 676.2663498, 1.0/Eps_rate*3.56,
        683.2639429, 1.0/Eps_rate*3.57, 690.331863, 1.0/Eps_rate*3.58, 697.4708169, 1.0/Eps_rate*3.59,
        704.6815185, 1.0/Eps_rate*3.6, 711.9646889, 1.0/Eps_rate*3.61, 719.3210563, 1.0/Eps_rate*3.62,
        726.7513564, 1.0/Eps_rate*3.63, 734.2563323, 1.0/Eps_rate*3.64, 741.8367345, 1.0/Eps_rate*3.65,
        749.493321, 1.0/Eps_rate*3.66, 757.2268574, 1.0/Eps_rate*3.67, 765.0381172, 1.0/Eps_rate*3.68,
        772.9278815, 1.0/Eps_rate*3.69, 780.8969391, 1.0/Eps_rate*3.7, 788.9460872, 1.0/Eps_rate*3.71,
        797.0761305, 1.0/Eps_rate*3.72, 805.2878822, 1.0/Eps_rate*3.73, 813.5821633, 1.0/Eps_rate*3.74,
        821.9598033, 1.0/Eps_rate*3.75, 830.42164, 1.0/Eps_rate*3.76, 838.9685196, 1.0/Eps_rate*3.77,
        847.6012967, 1.0/Eps_rate*3.78, 856.3208347, 1.0/Eps_rate*3.79, 865.1280055, 1.0/Eps_rate*3.8,
        874.0236899, 1.0/Eps_rate*3.81, 883.0087773, 1.0/Eps_rate*3.82, 892.0841664, 1.0/Eps_rate*3.83,
        901.2507647, 1.0/Eps_rate*3.84, 910.5094888, 1.0/Eps_rate*3.85, 919.8612646, 1.0/Eps_rate*3.86,
        929.3070274, 1.0/Eps_rate*3.87, 938.8477216, 1.0/Eps_rate*3.88, 948.4843014, 1.0/Eps_rate*3.89,
        958.2177305, 1.0/Eps_rate*3.9, 968.0489821, 1.0/Eps_rate*3.91, 977.9790395, 1.0/Eps_rate*3.92,
        988.0088956, 1.0/Eps_rate*3.93, 998.1395534, 1.0/Eps_rate*3.94, 1008.372026, 1.0/Eps_rate*3.95,
        1018.707337, 1.0/Eps_rate*3.96, 1029.146519, 1.0/Eps_rate*3.97, 1039.690617, 1.0/Eps_rate*3.98,
        1050.340685, 1.0/Eps_rate*3.99, 1061.097787, 1.0/Eps_rate*4.0, 1071.963001]

    Strain_Table = []
    [div,rest] = np.divmod(len(Strain_Time_Table),8)   
    for i in range(div):
        Strain_Table.append(Strain_Time_Table[i*8:i*8+8])     
    Strain_Table.append(Strain_Time_Table[div*8:div*8+rest])
    
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('*Amplitude, name=Amp-1\n')
        for line in Strain_Table:
            str_line = ', '.join([' ' + str(round(x, 7)) for x in line])
            inp_file.write(str_line + '\n')
            
  
def material(file_name, E, mu, rho, Plastic, law):
    
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('**\n')
        inp_file.write('** ---------------------------------------------------------------------\n')
        inp_file.write('** MATERIALS\n')
        inp_file.write('**\n')
        inp_file.write('*Material, name=Material-1\n')
        inp_file.write('*Density\n')
        inp_file.write(' %.3e\n'%rho)
        inp_file.write('**\n')
        inp_file.write('*Elastic\n')
        inp_file.write(' %.2f, %.4f\n'%(E, mu))
        inp_file.write('**\n')
        inp_file.write('*Plastic\n')
        for line in Plastic:
            str_line = ', '.join([' ' + str(round(x, 8)) for x in line])
            inp_file.write(str_line + '\n')
        inp_file.write('**\n')  
            
        law_params_list = []
        if YIELD_FUNC_LOOKUP.get(law['type']):
            law_params_list = YIELD_FUNC_LOOKUP.get(law['type'])(law['parameters'])        
        
        PARAMS_PER_LINE = 8
        law_params_list_str = ''
        for idx, i in enumerate(law_params_list):
            if (idx % PARAMS_PER_LINE) == (PARAMS_PER_LINE - 1):
                end_char = '\n'
            elif idx == (len(law_params_list) - 1):
                end_char = '\n'
            else:
                end_char = ', '
            law_params_list_str += f'{i:.6f}{end_char}'

        if law['type'] == 'Tresca':
            inp_file.write('*potential, type=tresca\n')
        elif law['type'] == 'Hill1948':
            inp_file.write('*potential, type=hill\n')
        elif law['type'] == 'Hosford':
            inp_file.write(f'*potential, type=hosford, power={law["power"]}\n')            
        elif law['type'] == 'Barlat_Yld91':
            inp_file.write(f'*potential, type=barlat91, power={law["power"]}\n')
        elif law['type'] == 'Barlat_Yld2004_18p':
            inp_file.write(f'*potential, type=barlat, power={law["power"]}\n')
        elif law['type'] != 'VonMises':
            raise ValueError(f'Unknown yield function type: {law["type"]}.')

        inp_file.write(law_params_list_str)

  
def step(file_name, time_step, dt_i, dt_min, dt_max, b1, b2, U_left, U_right,
         U_up, U_bottom, Max_plastic_strain, num_interval):
    
    ###################################
    ####           Step           #####
    ###################################    
    # Default Abaqus value if 0
    if dt_i == 0:
        dt_i = ''
        
    if dt_min == 0:
        dt_min = '' 
    
    if dt_max == 0:
        dt_max = ''
      
    with open(file_name, 'a') as inp_file:
        inp_file.write('**\n')
        inp_file.write('** ---------------------------------------------------------------------\n')
        inp_file.write('** STEP\n')
        inp_file.write('**\n')     
        inp_file.write('*Step, name=Step-1, nlgeom=YES\n')
        inp_file.write('*Dynamic, Explicit\n')
        inp_file.write('%s, %s, %s, %s\n'%(str(dt_i), str(time_step), str(dt_min), str(dt_max))) 
        inp_file.write('*Bulk Viscosity\n')
        inp_file.write(' %.3f, %.3f\n'%(b1, b2))
        inp_file.write('**\n')
        inp_file.write('**\n')

        ###################################
        ####           BC             #####
        ###################################  
        inp_file.write('** BOUNDARY CONDITIONS\n')
        inp_file.write('**\n') 
        
        if U_left == 'free' and U_right == 'free':
            U_left = 0.0
        if U_bottom == 'free' and U_up == 'free':
            U_bottom = 0.0       	   	
        if  U_left != 'free':
            inp_file.write('** Name: BC-1 (Left) Type: Displacement/Rotation\n')
            inp_file.write('*Boundary, amplitude=Amp-1\n')
            inp_file.write('Left_nodes, 1, 1, %.5f\n'%U_left)
        if U_right != 'free':
            inp_file.write('** Name: BC-1 (Right) Type: Displacement/Rotation\n')
            inp_file.write('*Boundary, amplitude=Amp-1\n')
            inp_file.write('Right_nodes, 1, 1, %.5f\n'%U_right)
        if  U_up != 'free':
            inp_file.write('** Name: BC-1 (Up) Type: Displacement/Rotation\n')
            inp_file.write('*Boundary, amplitude=Amp-1\n')
            inp_file.write('Upper_nodes, 2, 2, %.5f\n'%U_up)
        if  U_bottom != 'free':
            inp_file.write('** Name: BC-1 (Bottom) Type: Displacement/Rotation\n')
            inp_file.write('*Boundary, amplitude=Amp-1\n')
            inp_file.write('Bottom_nodes, 2, 2, %.5f\n'%U_bottom)
        inp_file.write('**\n')
        inp_file.write('**\n')
        
        ###################################
        ####    Stopping criteria     #####
        ################################### 
        if Max_plastic_strain > 0.0:  
	        inp_file.write('*EXTREME VALUE, HALT=YES\n')
	        inp_file.write('*EXTREME ELEMENT VALUE, ELSET=MK_sample-1.Sheet_Elements, MAX\n')    
	        inp_file.write('\n')
	        inp_file.write('PEEQ, %.4f\n'%Max_plastic_strain)        
	        inp_file.write('**\n')
	        inp_file.write('**\n')        
        ###################################
        ####           Outputs        #####
        ###################################  
        inp_file.write('** OUTPUT REQUESTS\n')
        inp_file.write('**\n')
        inp_file.write('*Restart, write, number interval=1, time marks=NO\n')
        inp_file.write('**\n')
        
        num_interval_str = f', NUMBER INTERVAL={num_interval}' if num_interval else ''
        
        inp_file.write(f'*Output, field{num_interval_str}\n')
        inp_file.write('*Node Output\n')
        inp_file.write('COORD, U, V\n')
        inp_file.write('*Element Output, directions=YES\n')
        inp_file.write('E, LE, PE, PEEQ, PEMAG, S, COORD\n')
        inp_file.write('**\n')
        inp_file.write('*Output, history, variable=PRESELECT\n')
        inp_file.write('**\n')
        inp_file.write('*End Step')

