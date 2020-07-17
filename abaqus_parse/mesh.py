import numpy as np
from more_itertools import pairwise
from vecmaths import geometry

def donut_polar_grid_quad(inner_r, outer_r, num_rays, num_arcs, quad=1,     include_borders=True):
    """
    Create a donut polar grid of nodes for one quadrant.
    
    Parameters
    ----------
    inner_r : float
        Size of inner radius of donut
    outer_r : float
        Size of outer radius of donut
    num_rays : int
        Number of rays within a quadrant of the circle (splits angles)
    num_arcs : int
        Number of arcs within a quadrant of the circle (splits radius)
    quad : int
        Circle quadrant to be populated, accepted values are 1, 2, 3 or 4 (default is 1)
    include_borders : bool
        Include border nodes. Default value is True.    
    
    Outputs
    -------
    tuple (ndarray, ndarray), shape=(N,)
        Polar coordinates distance-angle for nodes.
        
    
    """
    if quad == 1:
        ang_lim_low = 0.0
        ang_lim_upp = np.pi / 2 + 0.1
        
    r = np.linspace(inner_r, outer_r, num_arcs )
    theta = np.linspace(0, 2 * np.pi, (num_rays - 1) * 4 +1)

    radius_matrix, theta_matrix = np.meshgrid(r,theta)
    quat_idx = np.intersect1d(np.where(theta_matrix.flatten() >= ang_lim_low)[0],
               np.where(theta_matrix.flatten() <= ang_lim_upp)[0])
    
    
    return theta_matrix.flatten()[quat_idx], radius_matrix.flatten()[quat_idx]


def order_nodes_elem(elem_nodes_all, all_nodes, all_nodes_labs):
    """
    Order the nodes of an element in anticlockwise direction, consistent with the element definition in Abaqus.
    
    Parameters:
    elem_nodes_all : array of (M, n)
        The nodes that form part of M elements with n corners each.
    all_nodes : array of (N, 2)
        The coordinates of all N nodes in the mesh.
    all_nodes_labs : array of (N,)
        The labels of all nodes in the mesh.
    
    """
    elem_nodes_all_ordered = []
    for el_idx, elem_nodes in enumerate(elem_nodes_all):
        
        
        # # sort nodes in growing order
        # node_el_sorted = np.sort(elem_nodes)
        
        # indices of nodes in node list
        node_idx = np.where(np.isin(all_nodes_labs, elem_nodes))[0]
       
        if len(node_idx)==3:
            vals, counts = np.unique(elem_nodes, return_counts=True)
            rep_idx = np.where(counts == 2)[0]
            node_idx = np.insert(node_idx, rep_idx, node_idx[np.where(counts == 2)])

            
        node_el_sorted_by_node_list = all_nodes_labs[node_idx]    
        # coords of nodes in the order they are in the all nodes list?
        node_coord_sorted = all_nodes[node_idx]
        # order nodes in correct order
        
        node_ord = geometry.order_coplanar_points(
            np.hstack((node_coord_sorted, np.zeros((len(node_coord_sorted),1)))).T,
            np.array([[0,0,1]]).T
        )
        
        elem_nodes_all_ordered.append(node_el_sorted_by_node_list[node_ord])
        
    return np.concatenate(elem_nodes_all_ordered).reshape(-1, 4)



def cells_from_nodes(nodes_array):
    """
    Find the cell vertices in a regular grid defined by an array of nodes
    
    """
    c = 0
    cells = []

    for ri, row in enumerate(nodes_array):
        count = 0
        if ri != len(nodes_array)-1:
            cells.append([])
        for item1, item2 in pairwise(row):
            if ri == 0:
                cells[ri].append([item1, item2])
            elif ri == len(nodes_array)-1:
                cells[ri-1][count].append(item1)
                cells[ri-1][count].append(item2)
    
            else:
                cells[ri].append([item1, item2])
                cells[ri-1][count].append(item1)
                cells[ri-1][count].append(item2)
            count += 1
    
    return np.array(cells)


def make_donut_mesh(rin, elem_side = 0.025, an=3,
                     sbox_width_mult_rad = 5, ret_rin_nodes=False):
    """
    
    Parameters
    ----------
    rin : float
        Inner radius in micrometers
    
    an : integer
     number of arcs (TO CALCULATE LATER) = 3 
    sbox_width_mult_rad : integer
       aproximate width of fan mesh in multiples of r 
    ret_rin_nodes : bool
        Return the inner nodes compromising the inner radius in anti-clockwise direction. 
       
    Returns
    -------
    
    
    """    
    rad_mm = rin * 1e-3
    
    # width of fan mesh in multiples of elem side
    fbox_side_mult_elem = int((sbox_width_mult_rad*rad_mm) / elem_side)
    fbox_side_mult_elem = round(fbox_side_mult_elem/3)*3
    
    # side size of fan mesh in mm
    fbox_side_mm = fbox_side_mult_elem * elem_side
    
    # number of rays
    rn = 2 * fbox_side_mult_elem + 1

    thetaout, rout  = donut_polar_grid_quad(rad_mm, fbox_side_mm,
                                            rn, an)
    
    xout, yout = geometry.polar2cart_2D(rout, thetaout)
    
    # Move nodes at ray ends to the sides of a square
    ##################################################
    # coords from 0 to an of side square box in multiples of ctip radii
    # half coordinates of the end of fan rays crossing square
    end_ray_half_coords =  np.linspace(0, fbox_side_mm, int((rn-1)/2+1))
    # print('end_ray_half_coords: ', end_ray_half_coords)

    end_ray_repeat_coords = (fbox_side_mm * np.ones(rn))[None]
    end_ray_coords = np.concatenate((
        np.concatenate((end_ray_half_coords[:-1],
                    np.flip(end_ray_half_coords)))[None],
            end_ray_repeat_coords
        ), axis=0)

    # The second half of array indices: 
    # 1st half + mid element + next element - 1 for zero based
    # int((rn-1)/2) + 1 + 1 - 1 
    flip_idx = int((rn-1)/2) + 1 - 1
    end_ray_coords[:, :flip_idx] = np.flip(end_ray_coords[:, :flip_idx], axis=0)

    # Array of nodes as (Rows, Cols, (x, y)_node)
    don_reshaped = np.concatenate((xout[None].T, yout[None].T), axis=1).reshape(rn, an, 2)
    don_reshaped[:,-1,:] = end_ray_coords.T
    don_nodes_flattened = don_reshaped.reshape(rn*an, 2)
    don_node_labs = np.arange(1, an * rn + 1).reshape(rn, an)
    print('don_node_labs: ', don_node_labs[:,0])
    don_node_labs_flattened = don_node_labs.flatten()
    
    
    don_cells = cells_from_nodes(don_node_labs)
    don_cells_flattened = don_cells.reshape(-1, 4)
    

    don_cells_flattened = order_nodes_elem(don_cells_flattened, don_nodes_flattened,don_node_labs_flattened)

    don_cell_centres = np.mean(
        cells_from_nodes(don_reshaped), axis=2
    )
    don_cell_centres_flattened = don_cell_centres.reshape((an-1) * (rn-1), 2)

    don_cell_labs = np.arange(1, (an-1) * (rn-1) + 1).reshape(an-1, rn-1)
    don_cell_labs_flattened = don_cell_labs.flatten()

    if ret_rin_nodes:
        return (don_nodes_flattened, don_node_labs_flattened,
            don_cells_flattened, don_cell_centres_flattened, don_cell_labs_flattened, don_node_labs[:,0])    
    else:
        return (don_nodes_flattened, don_node_labs_flattened,
            don_cells_flattened, don_cell_centres_flattened, don_cell_labs_flattened)

def make_fine_plus_donut_mesh(rin, ref_mesh_side=0.2, elem_side = 0.025, an=4, sbox_width_mult_rad=5, ret_crack_line=True):
    """
    Build refined mesh = 4 rectilinear + 1 donut mesh at crack tip
    ''''''''''''''''''''''''''''''''''''''
    '               '    '               '
    '        5      '  4 '        3      '
    '               '    '               '
    '               '    '               '
    '               '    '               '
    '               '    '               '
    ''''''''''''''''''''''''''''''''''''''
    '        6      ' d-1'        2      '
    ''''''''''''''''' '  '               '
                      ''''''''''''''''''''
    
    
    """
    rad_mm = rin * 1e-3
    # number of elements on side of fine mesh
    num_elem_side = int(ref_mesh_side / elem_side)
    line_elem =  num_elem_side - 1 #  substract corner element
    line_elem = round(line_elem / 3) * 3
    num_elem_side = line_elem + 1
    ref_mesh_size = [2 * num_elem_side * elem_side, num_elem_side * elem_side]
    
    #refined mesh nodes and cells
    ref_mesh_nodes_all = []
    ref_mesh_labs_all = []
    ref_mesh_cells_all = []
    ref_mesh_cell_centres_all = []
    ref_mesh_cell_labs_all = []
    
    # donut mesh nodes and cells
    if ret_crack_line:
        don_nodes_flattened, don_node_labs_flattened,\
    don_cells, don_cell_centres_flattened, don_cell_labs_flattened, crack_line =\
    make_donut_mesh(rin, elem_side=elem_side, an=an, 
                    sbox_width_mult_rad=sbox_width_mult_rad, ret_rin_nodes=True)
    else:
        don_nodes_flattened, don_node_labs_flattened,\
    don_cells, don_cell_centres_flattened, don_cell_labs_flattened =\
    make_donut_mesh(rin, elem_side=elem_side, an=an, 
                    sbox_width_mult_rad=sbox_width_mult_rad)
        
    rn = int(don_nodes_flattened.shape[0] / an)
    don_mesh_size = [
           ( don_nodes_flattened.reshape(rn, an, 2))[0][-1][0],
            (don_nodes_flattened.reshape(rn, an, 2))[0][-1][0]
    ]
    
    ref_mesh_nodes_all.append(don_nodes_flattened)
    ref_mesh_labs_all.append(don_node_labs_flattened)
    ref_mesh_cells_all.append(don_cells)
    ref_mesh_cell_centres_all.append(don_cell_centres_flattened)
    ref_mesh_cell_labs_all.append(don_cell_labs_flattened)


    rect_mesh_shifts = [
        [don_mesh_size[0], 0],
        [don_mesh_size[0], don_mesh_size[1]],
        [0, don_mesh_size[1]],
        [-ref_mesh_size[0] / 2, don_mesh_size[1]],
        [-ref_mesh_size[0] / 2, 0],
    ]
    rect_mesh_sizes = [
        [ref_mesh_size[0] / 2 - don_mesh_size[0], don_mesh_size[1]],
        [ref_mesh_size[0] / 2 - don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] / 2, ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] / 2, don_mesh_size[1] - rad_mm]
    ]

    for i in range(5):
        mesh_size = rect_mesh_sizes[i]
        mesh_grid = [int(mesh_size[0]/elem_side) + 1,
                     int(mesh_size[1]/elem_side) + 1]

        w = np.linspace(0, mesh_size[0], mesh_grid[0])  
        h = np.linspace(0, mesh_size[1], mesh_grid[1])

        if i == 4:
            h = list(don_nodes_flattened.reshape(rn, an, 2)[-1, :, 1])  # DONT have this
            mesh_grid[1] = len(h)
        meshw, meshh = np.meshgrid(w, h)

        # shift meshgrid 
        meshw += rect_mesh_shifts[i][0]
        meshh += rect_mesh_shifts[i][1]

        # x, y coordinates of nodes
        nodes_coord = np.concatenate([meshw.flatten()[None],
                                      meshh.flatten()[None]], axis=0).T
        node_labs = np.zeros_like(meshw.flatten())
        reshaped = nodes_coord.reshape(mesh_grid[1], mesh_grid[0], 2)
        
        for ni, n in enumerate(nodes_coord):
            # find which nodes already been created 
            cond = np.isclose(n, np.concatenate(ref_mesh_nodes_all), 0.0001)
            found_idx = np.where(np.all(cond, axis=1))        
            if found_idx[0].size>0:
                node_labs[ni] = np.concatenate(ref_mesh_labs_all)[found_idx]

        # Apply mask over repeated border nodes
        mesh_mask = np.zeros_like(meshw)
        if i == 0:
            mesh_mask[:,0] = 1
        elif i == 1:
            mesh_mask[0,:] = 1
        elif i == 2:
            mesh_mask[0,:] = 1
            mesh_mask[:,-1] = 1
        elif i == 3:
            mesh_mask[:,-1] = 1
        elif i == 4:
            mesh_mask[-1,:] = 1
            mesh_mask[:,-1] = 1

        meshw, meshh = [np.ma.masked_array(i, mask=mesh_mask) for i in np.meshgrid(w, h)]
        # shift meshgrid 
        meshw += rect_mesh_shifts[i][0]
        meshh += rect_mesh_shifts[i][1]
        nodes_coord_masked = np.concatenate([meshw.compressed()[None],
                                      meshh.compressed()[None]], axis=0).T

        new_node_idx = np.where(node_labs==0)
        node_labs[new_node_idx] = np.arange(ref_mesh_labs_all[-1][-1]+1,
                                            ref_mesh_labs_all[-1][-1]+1 + len(meshw.compressed()))
        node_labs_grid = node_labs.reshape(meshw.shape)
        
        ref_mesh_nodes_all.append(nodes_coord_masked)
        ref_mesh_labs_all.append(node_labs[new_node_idx])

        # CHANGE!
        cells = (cells_from_nodes(node_labs_grid)).reshape(-1, 4)
        cells = order_nodes_elem(cells, np.concatenate(ref_mesh_nodes_all),
                                 np.concatenate(ref_mesh_labs_all))
        cell_centres = np.mean(
            cells_from_nodes(reshaped), axis=2
        )
        cell_centres_flattened = cell_centres.reshape((mesh_grid[0]-1) * (mesh_grid[1]-1), 2)
        cell_labs = np.arange(ref_mesh_cell_labs_all[-1][-1]+1,
                                            ref_mesh_cell_labs_all[-1][-1]+1 + (mesh_grid[0]-1) * (mesh_grid[1]-1))

        ref_mesh_cells_all.append(cells.astype('int'))
        ref_mesh_cell_centres_all.append(cell_centres_flattened)
        ref_mesh_cell_labs_all.append(cell_labs)
    if ret_crack_line:
        return ref_mesh_nodes_all, ref_mesh_labs_all,\
            ref_mesh_cells_all, ref_mesh_cell_centres_all, ref_mesh_cell_labs_all, crack_line
    else:    
        return ref_mesh_nodes_all, ref_mesh_labs_all,\
            ref_mesh_cells_all, ref_mesh_cell_centres_all, ref_mesh_cell_labs_all

def find_border(nodes, node_labels, condition):
    """
    Parameters
    ----------
    nodes: ndarray of (N, 2)
    condition: string
        One of minx, maxx, miny, maxy
    
    """
    if condition=='miny':
        nodes_bord = np.min(nodes[:,1])
        nodes_bord_idx =  np.where(nodes[:,1]==nodes_bord)[0]
    elif condition=='maxy':
        nodes_bord = np.max(nodes[:,1])
        nodes_bord_idx =  np.where(nodes[:,1]==nodes_bord)[0]
    elif condition=='minx':
        nodes_bord = np.min(nodes[:,0])
        nodes_bord_idx =  np.where(nodes[:,0]==nodes_bord)[0]
    elif condition=='maxx':
        nodes_bord = np.max(nodes[:,0])
        nodes_bord_idx =  np.where(nodes[:,0]==nodes_bord)[0]
    
    return nodes[nodes_bord_idx], node_labels[nodes_bord_idx]


def expand_mesh(coords, coord_labs, expand_dir, max_node_label, max_elem_label, exp_modes=['transition'],
                num_layers=[1], corner_beg=False, corner_end=False, end_dir=None, 
                width=1, width_mult=1.5):
    """
    Parameters
    ----------
    width : float
        Width of a single layer [length].
    
    """
    
    nodes_all = []
    nodes_all_labels = []
    elements_all = []
    elements_all_labels = []
    corner_node_labs = []
    # keep track of the outer most nodes and labels
    coords_outer = np.copy(coords)
    labs_outer = np.copy(coord_labs)
    width_outer = width
    max_node_label = max_node_label
    
    for mi, mode in enumerate(exp_modes):
        n = len(coords_outer)
        vec_outer = expand_dir * width_outer * num_layers[mi]
        if mode=='transition':
            
            # number of blocks of three-to-one elements = 
            nb = (n - 1) // 3
            # number of new outside nodes = N
            N = nb + 1 
            M = (n - 1) % 3      # number of extra nodes (0, 1 or 2)

            vec_inner = 0.5 * vec_outer

            coords_block = coords_outer[:3*nb+1]
            coord_labs_block = labs_outer[:3*nb+1]

            # Expand main line 
            nodes_outer = coords_block[::3] + vec_outer
            nodes_outer_labels = np.arange(max_node_label+1, 
                                            max_node_label+len(nodes_outer)+1)
            max_node_label += len(nodes_outer)
            nodes_inner_a = coords_block[1::3] + vec_inner
            nodes_inner_a_labels = np.arange(max_node_label+1, max_node_label+len(nodes_inner_a)+1)
            max_node_label += len(nodes_inner_a)
            nodes_inner_b = coords_block[2::3] + vec_inner
            nodes_inner_b_labels = np.arange(max_node_label+1, max_node_label+len(nodes_inner_b)+1)
            max_node_label += len(nodes_inner_b)
            
            
            nodes_all_layer = np.concatenate((
                nodes_outer,
                nodes_inner_a,
                nodes_inner_b,
            ))
            nodes_all_labels_layer = np.concatenate((
                nodes_outer_labels,
                nodes_inner_a_labels,
                nodes_inner_b_labels,
            ))
            assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
            
            # find elements
            elem_A = np.array([
                coord_labs_block[:-1:3], coord_labs_block[1::3], 
                nodes_inner_a_labels, nodes_outer_labels[:-1]]).T

            elem_B = np.array([
                coord_labs_block[1::3], coord_labs_block[2::3],
                nodes_inner_a_labels, nodes_inner_b_labels
            ]).T

            elem_C = np.array([
                coord_labs_block[2::3], coord_labs_block[::3][1:],
                nodes_inner_b_labels, nodes_outer_labels[1:]
            ]).T

            elem_D = np.array([
                nodes_outer_labels[:-1], nodes_outer_labels[1:],
                nodes_inner_a_labels, nodes_inner_b_labels
            ]).T
            
            elements_all_layer = np.concatenate([
                elem_A,
                elem_B,
                elem_C,
                elem_D
            ])
   
            if M > 0 or corner_end:
                node_end = coords_outer[-1] + vec_outer
                
                
                node_end_lab = max_node_label + 1
                if corner_end:
                    node_end += end_dir * width_outer
                    corner_node_labs.append(node_end_lab)
                    
                # add node to all nodes in layer
                nodes_all_layer = np.concatenate((nodes_all_layer, node_end[None]))
                
                nodes_all_labels_layer = np.concatenate((nodes_all_labels_layer, [node_end_lab]))
                assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
                
                # add node to outer layer
                nodes_outer = np.concatenate((nodes_outer, node_end[None]))
                nodes_outer_labels = np.concatenate((nodes_outer_labels, [node_end_lab]))
                
                if M != 2:        
                    elem_end = np.array([labs_outer[-1], coord_labs_block[-1], 
                                node_end_lab, nodes_outer_labels[-2]])[None]
                    elements_all_layer = np.append(elements_all_layer, elem_end, axis=0)
                
                elif M == 2:
                    node_extra = np.mean((nodes_outer[-2], node_end), axis=0)
                    node_extra_lab = node_end_lab + 1
                    nodes_all_layer = np.concatenate((nodes_all_layer, node_extra[None]))
                    
                    nodes_all_labels_layer = np.concatenate((nodes_all_labels_layer, [node_extra_lab]))
                    assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
                    
                    nodes_outer = np.concatenate((nodes_outer, node_extra[None]))
                    nodes_outer_labels = np.concatenate((nodes_outer_labels, [node_extra_lab]))

                    elem_extra = np.array([
                        [labs_outer[-2],  labs_outer[-3],
                            nodes_outer_labels[-3], node_extra_lab],
                        [labs_outer[-2], labs_outer[-1],
                            node_extra_lab, node_end_lab]
                    ])
                    elements_all_layer = np.append(elements_all_layer, elem_extra, axis=0)
                    
            if corner_beg:
                node_beg = coords_outer[0] + vec_outer - end_dir * width_outer
                nodes_all_layer = np.concatenate((nodes_all_layer, node_beg[None]))
                node_beg_label = nodes_all_labels_layer[-1]+1
                
                nodes_all_labels_layer = np.append(nodes_all_labels_layer, node_beg_label)
                assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
                
                nodes_outer = np.concatenate((nodes_outer, node_beg[None]))
                nodes_outer_labels = np.concatenate((nodes_outer_labels, [node_beg_label]))
                
                elem_beg = np.array([node_beg_label, coord_labs_block[0], 
                                coord_labs_block[0], nodes_outer_labels[0]])[None]
                elements_all_layer = np.append(elements_all_layer, elem_beg, axis=0)
            

        elif mode=='uniform':
            nodes_outer = coords_outer + vec_outer
            nodes_outer_labels = np.arange(max_node_label+1, 
                                            max_node_label+len(nodes_outer)+1)
            if corner_beg:
                node_beg = coords_outer[0] + vec_outer - end_dir * width_outer * num_layers[mi]
                nodes_outer[0] = node_beg
            if corner_end:
                nodes_outer[-1] = nodes_outer[-1] + end_dir * width_outer * num_layers[mi]
                
            nodes_all_layer = np.copy(nodes_outer)
            
            nodes_all_labels_layer = np.copy(nodes_outer_labels)
            assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
            
            layers_inner = np.arange(1, num_layers[mi], 1)
        
            max_node_lab_inner = nodes_outer_labels[-1]
            if layers_inner.size > 0:
                nodes_inner = []
                nodes_inner_labels = []
                for layer in layers_inner:
                    vec_inner = expand_dir * width_outer * layer
                    nodes_inner_i = coords_outer + vec_inner
                    if corner_beg:
                        nodes_inner_i[0] = nodes_inner_i[0] - end_dir * width_outer * layer
                    if corner_end:
                        nodes_inner_i[-1] = nodes_inner_i[-1] + end_dir * width_outer * layer
                    nodes_inner.append(nodes_inner_i)
                    nodes_inner_labels.append(
                        np.arange(
                            max_node_lab_inner+1, 
                            max_node_lab_inner+len(coords_outer)+1
                        )
                    )

                    max_node_lab_inner = np.max(np.array(nodes_inner_labels))

                nodes_all_layer = np.concatenate((nodes_all_layer, np.squeeze(np.concatenate(nodes_inner))))
                
                nodes_all_labels_layer = np.append(nodes_all_labels_layer, nodes_inner_labels)
                assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)
            
            nodes_all_labels_layer_reshaped = nodes_all_labels_layer.reshape(len(nodes_all_labels_layer)//len(coords_outer), len(coords_outer))
            
            labs_outer_sort = labs_outer[coords_outer[:,np.where(end_dir==1)[0][0]].argsort()]
            elem_nodes_concat = np.concatenate((nodes_all_labels_layer_reshaped, labs_outer_sort[None]))
            elements_all_layer = (cells_from_nodes(elem_nodes_concat)).reshape(-1, 4)
            

        elements_all_labels_layer = np.arange(max_elem_label+1, 
                                            max_elem_label+len(elements_all_layer)+1)

        coords_outer = np.copy(nodes_outer)
        labs_outer = np.copy(nodes_outer_labels)
        labs_outer = labs_outer[coords_outer[:,np.where(end_dir==1)[0][0]].argsort()]
        coords_outer = coords_outer[coords_outer[:,np.where(end_dir==1)[0][0]].argsort()]
            
        nodes_all.append(nodes_all_layer)
        nodes_all_labels.append(nodes_all_labels_layer)
        elements_all.append(elements_all_layer)
        elements_all_labels.append(elements_all_labels_layer)
        
        width_outer = (mi + 1) * width_mult * width
        max_node_label = np.max(np.concatenate(nodes_all_labels))
        max_elem_label = np.max(elements_all_labels_layer)

    
    elements_all = order_nodes_elem(np.concatenate(elements_all).astype('int'),
                                    np.concatenate((coords, np.concatenate(nodes_all))), 
                                 np.concatenate((coord_labs, np.concatenate(nodes_all_labels))))
    return np.concatenate(nodes_all), np.concatenate(nodes_all_labels),\
            elements_all, np.concatenate(elements_all_labels)

def expand_mult_dirs(conds, expand_dirs, end_dirs, mesh_nodes_all, 
                    mesh_node_labs_all, mesh_elem_all, mesh_elem_labs_all,
                    num_layers=[3,3], width=1):
    
    for i in range(len(conds)):

        # find max node and element labels to start
        max_node_label = np.max(mesh_node_labs_all)
        max_elem_label = np.max(mesh_elem_labs_all)


        bord_nodes, bord_labs = find_border(
            mesh_nodes_all, 
            mesh_node_labs_all,
            conds[i]
        )    
        sort_idx_x = np.lexsort((bord_nodes[:,1], bord_nodes[:,0],))
        expand_dir = expand_dirs[i]
        end_dir = end_dirs[i]

        expmesh_nodes, expmesh_labels, elements_all, elements_all_labels = expand_mesh(
            bord_nodes[sort_idx_x], 
            bord_labs[sort_idx_x],
            expand_dir, max_node_label=max_node_label,
            max_elem_label=max_elem_label, 
            end_dir=end_dir, width=width,
            exp_modes=[ 'uniform',],
            num_layers=num_layers, 
        )

        mesh_nodes_all = np.concatenate(
            (mesh_nodes_all,
            expmesh_nodes)
        )
        mesh_node_labs_all = np.concatenate(
            (mesh_node_labs_all,
            expmesh_labels)
        )
        mesh_elem_all = np.concatenate(
            (mesh_elem_all,
            elements_all)
        )
        mesh_elem_labs_all = np.concatenate(
            (mesh_elem_labs_all,
            elements_all_labels)
        )
        
    return mesh_nodes_all, mesh_node_labs_all, mesh_elem_all.astype('int'), mesh_elem_labs_all.astype('int')