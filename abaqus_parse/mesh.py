import numpy as np
from more_itertools import pairwise
from abaqus_parse import fracture, utils


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

    r = np.linspace(inner_r, outer_r, num_arcs)
    theta = np.linspace(0, 2 * np.pi, (num_rays - 1) * 4 + 1)

    radius_matrix, theta_matrix = np.meshgrid(r, theta)
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

        if len(node_idx) == 3:
            vals, counts = np.unique(elem_nodes, return_counts=True)
            rep_idx = np.where(counts == 2)[0]
            node_idx = np.insert(node_idx, rep_idx, node_idx[np.where(counts == 2)])

        node_el_sorted_by_node_list = all_nodes_labs[node_idx]
        # coords of nodes in the order they are in the all nodes list?
        node_coord_sorted = all_nodes[node_idx]
        # order nodes in correct order

        node_ord = utils.order_coplanar_points(
            np.hstack((node_coord_sorted, np.zeros((len(node_coord_sorted), 1)))).T,
            np.array([[0, 0, 1]]).T
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


def make_donut_mesh(crack_tip_radius_microns, fine_mesh_element_length=0.025, fan_box_num_side_elements=3, fan_box_width=5):
    """

    Parameters
    ----------
    crack_tip_radius_microns : float
        Inner radius in micrometers

    fan_box_num_side_elements : integer
     number of arcs (TO CALCULATE LATER) = 3 
    fan_box_width : integer
       aproximate width of fan mesh in multiples of r 
    ret_rin_nodes : bool
        Return the inner nodes compromising the inner radius in anti-clockwise direction. 

    Returns
    -------


    """
    rad_mm = crack_tip_radius_microns * 1e-3

    # width of fan mesh in multiples of elem side
    fbox_side_mult_elem = int((fan_box_width*rad_mm) / fine_mesh_element_length)
    fbox_side_mult_elem = round(fbox_side_mult_elem/3)*3

    # side size of fan mesh in mm
    fbox_side_mm = fbox_side_mult_elem * fine_mesh_element_length

    # number of rays
    rn = 2 * fbox_side_mult_elem + 1

    thetaout, rout = donut_polar_grid_quad(rad_mm, fbox_side_mm,
                                           rn, fan_box_num_side_elements)

    xout, yout = utils.polar2cart_2D(rout, thetaout)

    # Move nodes at ray ends to the sides of a square
    ##################################################
    # coords from 0 to fan_box_num_side_elements of side square box in multiples of ctip radii
    # half coordinates of the end of fan rays crossing square
    end_ray_half_coords = np.linspace(0, fbox_side_mm, int((rn-1)/2+1))
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
    don_reshaped = np.concatenate((xout[None].T, yout[None].T), axis=1).reshape(
        rn, fan_box_num_side_elements, 2)
    don_reshaped[:, -1, :] = end_ray_coords.T
    don_nodes_flattened = don_reshaped.reshape(rn*fan_box_num_side_elements, 2)
    don_node_labs = np.arange(1, fan_box_num_side_elements *
                              rn + 1).reshape(rn, fan_box_num_side_elements)
    don_node_labs_flattened = don_node_labs.flatten()

    don_cells = cells_from_nodes(don_node_labs)
    don_cells_flattened = don_cells.reshape(-1, 4)

    don_cells_flattened = order_nodes_elem(
        don_cells_flattened, don_nodes_flattened, don_node_labs_flattened)

    don_cell_centres = np.mean(
        cells_from_nodes(don_reshaped), axis=2
    )
    don_cell_centres_flattened = don_cell_centres.reshape(
        (fan_box_num_side_elements-1) * (rn-1), 2)

    don_cell_labs = np.arange(1, (fan_box_num_side_elements-1)
                              * (rn-1) + 1).reshape(fan_box_num_side_elements-1, rn-1)
    don_cell_labs_flattened = don_cell_labs.flatten()

    out = {
        'nodes_flattened': don_nodes_flattened,
        'node_labels_flattened': don_node_labs_flattened,
        'cells_flattened': don_cells_flattened,
        'cell_centres_flattened': don_cell_centres_flattened,
        'cell_labels_flattened': don_cell_labs_flattened,
    }

    return out


def make_fine_plus_donut_mesh(crack_tip_radius_microns, fine_mesh_length=0.2, fine_mesh_element_length=0.025, fan_box_num_side_elements=4,  fan_box_width=5, ret_crack_definition=True, size_behind_crack=0.25):
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

    Parameters
    ----------
    size_behind_crack : float
        Size of the mesh behind the crack as fraction of fine_mesh_length.

    """
    rad_mm = crack_tip_radius_microns * 1e-3
    # number of elements on side of fine mesh
    num_elem_side = int(fine_mesh_length / fine_mesh_element_length)
    line_elem = num_elem_side - 1  # substract corner element
    line_elem = round(line_elem / 3) * 3
    num_elem_side = line_elem + 1
    ref_mesh_size = [2 * num_elem_side * fine_mesh_element_length,
                     num_elem_side * fine_mesh_element_length]

    # refined mesh nodes and cells
    ref_mesh_nodes_all = []
    ref_mesh_labs_all = []
    ref_mesh_cells_all = []
    ref_mesh_cell_centres_all = []
    ref_mesh_cell_labs_all = []

    # donut mesh nodes and cells
    donut_mesh = make_donut_mesh(
        crack_tip_radius_microns,
        fine_mesh_element_length=fine_mesh_element_length, fan_box_num_side_elements=fan_box_num_side_elements,
        fan_box_width=fan_box_width,
    )

    rn = int(donut_mesh['nodes_flattened'].shape[0] / fan_box_num_side_elements)
    don_mesh_size = [
        (donut_mesh['nodes_flattened'].reshape(
            rn, fan_box_num_side_elements, 2))[0][-1][0],
        (donut_mesh['nodes_flattened'].reshape(
            rn, fan_box_num_side_elements, 2))[0][-1][0]
    ]
    if ret_crack_definition:
        crack_line = donut_mesh['node_labels_flattened'].reshape(
            rn, fan_box_num_side_elements)[:, 0]
        crack_front = donut_mesh['node_labels_flattened'].reshape(
            rn, fan_box_num_side_elements)[0]

    ref_mesh_nodes_all.append(donut_mesh['nodes_flattened'])
    ref_mesh_labs_all.append(donut_mesh['node_labels_flattened'])
    ref_mesh_cells_all.append(donut_mesh['cells_flattened'])
    ref_mesh_cell_centres_all.append(donut_mesh['cell_centres_flattened'])
    ref_mesh_cell_labs_all.append(donut_mesh['cell_labels_flattened'])

    rect_mesh_shifts = [
        [don_mesh_size[0], 0],
        [don_mesh_size[0], don_mesh_size[1]],
        [0, don_mesh_size[1]],
        [-ref_mesh_size[0] * size_behind_crack, don_mesh_size[1]],
        [-ref_mesh_size[0] * size_behind_crack, 0],
    ]
    rect_mesh_sizes = [
        [ref_mesh_size[0] / 2 - don_mesh_size[0], don_mesh_size[1]],
        [ref_mesh_size[0] / 2 - don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] * size_behind_crack, ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] * size_behind_crack, don_mesh_size[1] - rad_mm]
    ]

    for i in range(5):
        mesh_size = rect_mesh_sizes[i]
        mesh_grid = [int(mesh_size[0]/fine_mesh_element_length) + 1,
                     int(mesh_size[1]/fine_mesh_element_length) + 1]

        w = np.linspace(0, mesh_size[0], mesh_grid[0])
        h = np.linspace(0, mesh_size[1], mesh_grid[1])

        if i == 4:
            h = list(donut_mesh['nodes_flattened'].reshape(
                rn, fan_box_num_side_elements, 2)[-1, :, 1])  # DON'T have this
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
            if found_idx[0].size > 0:
                node_labs[ni] = np.concatenate(ref_mesh_labs_all)[found_idx]

        # Apply mask over repeated border nodes
        mesh_mask = np.zeros_like(meshw)
        if i == 0:
            mesh_mask[:, 0] = 1
        elif i == 1:
            mesh_mask[0, :] = 1
        elif i == 2:
            mesh_mask[0, :] = 1
            mesh_mask[:, -1] = 1
        elif i == 3:
            mesh_mask[:, -1] = 1
        elif i == 4:
            mesh_mask[-1, :] = 1
            mesh_mask[:, -1] = 1

        meshw, meshh = [np.ma.masked_array(i, mask=mesh_mask) for i in np.meshgrid(w, h)]
        # shift meshgrid
        meshw += rect_mesh_shifts[i][0]
        meshh += rect_mesh_shifts[i][1]
        nodes_coord_masked = np.concatenate([meshw.compressed()[None],
                                             meshh.compressed()[None]], axis=0).T

        new_node_idx = np.where(node_labs == 0)
        node_labs[new_node_idx] = np.arange(ref_mesh_labs_all[-1][-1]+1,
                                            ref_mesh_labs_all[-1][-1]+1 + len(meshw.compressed()))
        node_labs_grid = node_labs.reshape(meshw.shape)
        if i == 0:
            # print('node_labs_grid[0,1:]: ', node_labs_grid[0,1:])
            crack_front = np.concatenate((crack_front, node_labs_grid[0, 1:]))
            # print('front: ', front)
        elif i == 4:
            crack_lip = node_labs_grid[0, :-1]
            # print('crack_lip: ', crack_lip)
        ref_mesh_nodes_all.append(nodes_coord_masked)
        ref_mesh_labs_all.append(node_labs[new_node_idx])

        # CHANGE!
        cells = (cells_from_nodes(node_labs_grid)).reshape(-1, 4)
        cells = order_nodes_elem(cells, np.concatenate(ref_mesh_nodes_all),
                                 np.concatenate(ref_mesh_labs_all))
        cell_centres = np.mean(
            cells_from_nodes(reshaped), axis=2
        )
        cell_centres_flattened = cell_centres.reshape(
            (mesh_grid[0]-1) * (mesh_grid[1]-1), 2)
        cell_labs = np.arange(ref_mesh_cell_labs_all[-1][-1]+1,
                              ref_mesh_cell_labs_all[-1][-1]+1 + (mesh_grid[0]-1) * (mesh_grid[1]-1))

        ref_mesh_cells_all.append(cells.astype('int'))
        ref_mesh_cell_centres_all.append(cell_centres_flattened)
        ref_mesh_cell_labs_all.append(cell_labs)
    crack_nodes = {
        'crack_lip': crack_lip.astype('int'),
        'crack_line': crack_line.astype('int'),
        'crack_front': crack_front.astype('int'),
    }

    out = {
        'node_coordinates': ref_mesh_nodes_all,
        'node_labels': ref_mesh_labs_all,
        'element_nodes': ref_mesh_cells_all,
        'element_centre_coordinates': ref_mesh_cell_centres_all,
        'element_labels': ref_mesh_cell_labs_all,
    }
    if ret_crack_definition:
        out['crack_definition'] = crack_nodes

    return out


def find_border(nodes, node_labels, condition):
    """
    Parameters
    ----------
    nodes: ndarray of (N, 2)
    condition: string
        One of minx, maxx, miny, maxy

    """
    if condition == 'miny':
        nodes_bord = np.min(nodes[:, 1])
        nodes_bord_idx = np.where(nodes[:, 1] == nodes_bord)[0]
    elif condition == 'maxy':
        nodes_bord = np.max(nodes[:, 1])
        nodes_bord_idx = np.where(nodes[:, 1] == nodes_bord)[0]
    elif condition == 'minx':
        nodes_bord = np.min(nodes[:, 0])
        nodes_bord_idx = np.where(nodes[:, 0] == nodes_bord)[0]
    elif condition == 'maxx':
        nodes_bord = np.max(nodes[:, 0])
        nodes_bord_idx = np.where(nodes[:, 0] == nodes_bord)[0]

    return nodes[nodes_bord_idx], node_labels[nodes_bord_idx]


def expand_mesh(coords, coord_labs, expand_dir, max_node_label, max_elem_label, exp_modes=['transition'], num_layers=[1], corner_beg=False, corner_end=False, end_dir=None, width=1, width_mult=1.5):
    """
    Parameters
    ----------
    coords : (N, 2)
    coord_labs : (N,)

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
        if mode == 'transition':

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
            nodes_inner_a_labels = np.arange(
                max_node_label+1, max_node_label+len(nodes_inner_a)+1)
            max_node_label += len(nodes_inner_a)
            nodes_inner_b = coords_block[2::3] + vec_inner
            nodes_inner_b_labels = np.arange(
                max_node_label+1, max_node_label+len(nodes_inner_b)+1)
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

                nodes_all_labels_layer = np.concatenate(
                    (nodes_all_labels_layer, [node_end_lab]))
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

                    nodes_all_labels_layer = np.concatenate(
                        (nodes_all_labels_layer, [node_extra_lab]))
                    assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)

                    nodes_outer = np.concatenate((nodes_outer, node_extra[None]))
                    nodes_outer_labels = np.concatenate(
                        (nodes_outer_labels, [node_extra_lab]))

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
                nodes_outer_labels = np.concatenate(
                    (nodes_outer_labels, [node_beg_label]))

                elem_beg = np.array([node_beg_label, coord_labs_block[0],
                                     coord_labs_block[0], nodes_outer_labels[0]])[None]
                elements_all_layer = np.append(elements_all_layer, elem_beg, axis=0)

        elif mode == 'uniform':
            nodes_outer = coords_outer + vec_outer
            nodes_outer_labels = np.arange(max_node_label+1,
                                           max_node_label+len(nodes_outer)+1)
            if corner_beg:
                node_beg = coords_outer[0] + vec_outer - \
                    end_dir * width_outer * num_layers[mi]
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
                        nodes_inner_i[0] = nodes_inner_i[0] - \
                            end_dir * width_outer * layer
                    if corner_end:
                        nodes_inner_i[-1] = nodes_inner_i[-1] + \
                            end_dir * width_outer * layer
                    nodes_inner.append(nodes_inner_i)
                    nodes_inner_labels.append(
                        np.arange(
                            max_node_lab_inner+1,
                            max_node_lab_inner+len(coords_outer)+1
                        )
                    )

                    max_node_lab_inner = np.max(np.array(nodes_inner_labels))

                nodes_all_layer = np.concatenate(
                    (nodes_all_layer, np.squeeze(np.concatenate(nodes_inner))))

                nodes_all_labels_layer = np.append(
                    nodes_all_labels_layer, nodes_inner_labels)
                assert len(set(nodes_all_labels_layer)) == len(nodes_all_labels_layer)

            nodes_all_labels_layer_reshaped = nodes_all_labels_layer.reshape(
                len(nodes_all_labels_layer)//len(coords_outer), len(coords_outer))

            labs_outer_sort = labs_outer[coords_outer[:,
                                                      np.where(end_dir == 1)[0][0]].argsort()]
            elem_nodes_concat = np.concatenate(
                (nodes_all_labels_layer_reshaped, labs_outer_sort[None]))
            elements_all_layer = (cells_from_nodes(elem_nodes_concat)).reshape(-1, 4)

        elements_all_labels_layer = np.arange(max_elem_label+1,
                                              max_elem_label+len(elements_all_layer)+1)

        coords_outer = np.copy(nodes_outer)
        labs_outer = np.copy(nodes_outer_labels)
        labs_outer = labs_outer[coords_outer[:, np.where(end_dir == 1)[0][0]].argsort()]
        coords_outer = coords_outer[coords_outer[:,
                                                 np.where(end_dir == 1)[0][0]].argsort()]

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
                     num_layers=[3, 3], width=1):

    for i in range(len(conds)):

        # find max node and element labels to start
        max_node_label = np.max(mesh_node_labs_all)
        max_elem_label = np.max(mesh_elem_labs_all)

        bord_nodes, bord_labs = find_border(
            mesh_nodes_all,
            mesh_node_labs_all,
            conds[i]
        )
        sort_idx_x = np.lexsort((bord_nodes[:, 1], bord_nodes[:, 0],))
        expand_dir = expand_dirs[i]
        end_dir = end_dirs[i]

        expmesh_nodes, expmesh_labels, elements_all, elements_all_labels = expand_mesh(
            bord_nodes[sort_idx_x],
            bord_labs[sort_idx_x],
            expand_dir, max_node_label=max_node_label,
            max_elem_label=max_elem_label,
            end_dir=end_dir, width=width,
            exp_modes=['uniform', ],
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


def transition_mesh(refined_mesh_definition):
    """
    Expand a refined mesh in required directions through a transition mesh with increasing element size. The following joining mesh feature is used:
     _____
    |_|_|_|
    | |_| |
    |/___\|

    Parameters
    ----------



    """

    ref_mesh_nodes_all = refined_mesh_definition['node_coordinates']
    ref_mesh_labs_all = refined_mesh_definition['node_labels']
    ref_mesh_cells_all = refined_mesh_definition['element_nodes']
    ref_mesh_cell_centres_all = refined_mesh_definition['element_centre_coordinates']
    ref_mesh_cell_labs_all = refined_mesh_definition['element_labels']
    crack_nodes = refined_mesh_definition['crack_definition']

    ref_mesh_nodes_all_con = np.concatenate(ref_mesh_nodes_all)

    conds = ['minx', 'maxy', 'maxx', ]

    min_x = np.min(ref_mesh_nodes_all_con[:, 0])
    max_x = np.max(ref_mesh_nodes_all_con[:, 0])
    max_y = np.max(ref_mesh_nodes_all_con[:, 1])

    # border idx of expansion directions
    ref_mesh_bord_idx = [
        np.where(ref_mesh_nodes_all_con[:, 0] == min_x)[0],
        np.where(ref_mesh_nodes_all_con[:, 1] == max_y)[0],
        np.where(ref_mesh_nodes_all_con[:, 0] == max_x)[0],
    ]

    # expand directions (normal to border)
    expand_dirs = [
        np.array([-1, 0]),
        np.array([0, 1]),
        np.array([1, 0])
    ]

    # start-end direction (parallel to border)
    end_dirs = [
        np.array([0, 1]),
        np.array([1, 0]),
        np.array([0, 1])
    ]
    # corner end beginning and end modes
    corner_end_modes = [True, True, True]
    corner_beg_modes = [False, True, False]

    # create copies of all mesh components as a start
    mesh_nodes_all = ref_mesh_nodes_all.copy()
    mesh_node_labs_all = ref_mesh_labs_all.copy()
    mesh_elem_all = ref_mesh_cells_all.copy()
    mesh_elem_labs_all = ref_mesh_cell_labs_all.copy()

    # find max node and element labels to start
    max_node_label = int(ref_mesh_labs_all[-1][-1])
    max_elem_label = int(ref_mesh_cell_labs_all[-1][-1])

    for i in range(len(conds)):
        bord_nodes, bord_labs = find_border(ref_mesh_nodes_all_con,
                                            np.concatenate(ref_mesh_labs_all), conds[i])
        sort_idx_x = np.lexsort((bord_nodes[:, 1], bord_nodes[:, 0],))
        expand_dir = expand_dirs[i]
        end_dir = end_dirs[i]

        expmesh_nodes, expmesh_labels, elements_all, elements_all_labels = expand_mesh(bord_nodes[sort_idx_x],
                                                                                       bord_labs[sort_idx_x],
                                                                                       expand_dir, max_node_label=max_node_label,
                                                                                       max_elem_label=max_elem_label,
                                                                                       corner_end=corner_end_modes[
                                                                                           i], corner_beg=corner_beg_modes[i],
                                                                                       end_dir=end_dir, width=0.1,
                                                                                       exp_modes=['transition', 'uniform', 'transition', 'uniform',
                                                                                                  'transition', 'uniform', ],
                                                                                       num_layers=[
                                                                                           1, 2, 1, 1, 1, 1],
                                                                                       )

        max_node_label = expmesh_labels[-1]
        max_elem_label = elements_all_labels[-1]
        mesh_nodes_all.append(
            expmesh_nodes
        )
        mesh_node_labs_all.append(
            expmesh_labels
        )
        mesh_elem_all.append(elements_all)
        mesh_elem_labs_all.append(elements_all_labels)

    mesh_nodes_all = np.concatenate(mesh_nodes_all)
    mesh_node_labs_all = np.concatenate(mesh_node_labs_all).astype('int')
    mesh_elem_all = np.concatenate(mesh_elem_all)
    mesh_elem_labs_all = np.concatenate(mesh_elem_labs_all)

    # REMOVE duplicates
    # find unique node coordinates
    mesh_nodes_all_uniq, index, counts = np.unique(
        mesh_nodes_all, axis=0,
        return_index=True,
        return_counts=True
    )

    # labels of unique nodes
    mesh_node_labs_all_uniq = mesh_node_labs_all[index]

    # indices of removed nodes
    rem_nodes_index = np.setdiff1d(np.arange(len(mesh_nodes_all)), index)
    removed_nodes = mesh_nodes_all[rem_nodes_index]

    mesh_elem_all_uniq = mesh_elem_all.copy()
    if np.any(counts > 2):
        raise Warning('Label repeated more than twice.')
    else:
        for ni, rep_node in enumerate(mesh_nodes_all_uniq[counts == 2]):

            pl_idx_ni = np.nonzero((mesh_nodes_all_uniq == rep_node).all(axis=1))
            rep_nodes_idx = np.nonzero((removed_nodes == rep_node).all(axis=1))
            rep_nodes = removed_nodes[rep_nodes_idx][0]
            rep_nodes_lab = mesh_node_labs_all[rem_nodes_index][rep_nodes_idx][0]
            mesh_elem_all_uniq[np.where(mesh_elem_all == rep_nodes_lab)
                               ] = mesh_node_labs_all_uniq[pl_idx_ni]

    out = {
        'node_coordinates': mesh_nodes_all_uniq,
        'node_labels': mesh_node_labs_all_uniq,
        'element_nodes': mesh_elem_all_uniq,
        'element_labels': mesh_elem_labs_all,
    }

    return out


def compact_tension_specimen_mesh(refined_mesh_definition, dimensions, size_type, fraction):

    transition_mesh_definition = transition_mesh(refined_mesh_definition)
    crack_nodes = refined_mesh_definition['crack_definition']

    nd = transition_mesh_definition['node_coordinates']
    ndlab = transition_mesh_definition['node_labels']
    el = transition_mesh_definition['element_nodes']
    ellab = transition_mesh_definition['element_labels']

    specimen = fracture.standard_specimen(
        size_type,
        dimensions=dimensions,
        fraction=fraction
    )

    tot_size_ref = (
        np.max(nd[:, 0]) - np.min(nd[:, 0]),
        np.max(nd[:, 1]) - np.min(nd[:, 1]),
    )

    exp_width = (specimen['A'] - tot_size_ref[0]) / 2
    exp_height = (specimen['E'] - tot_size_ref[1])

    num_lays_width = int(exp_width // 1)
    width_lays_width = 1 + exp_width % 1 / num_lays_width

    num_lays_height = int(exp_height // 1)
    width_lays_height = 1 + exp_height % 1 / num_lays_height

    for i in range(num_lays_width):
        # find idx of border nodes
        conds = ['minx', 'maxx']

        # expand directions (normal to border)
        expand_dirs = [
            np.array([-1, 0]),
            np.array([1, 0])
        ]

        # start-end direction (parallel to border)
        end_dirs = [
            np.array([0, 1]),
            np.array([0, 1])
        ]

        nd, ndlab, el, ellab = expand_mult_dirs(conds, expand_dirs, end_dirs,
                                                nd, ndlab, el, ellab,
                                                num_layers=[1, 1], width=width_lays_width)

    for i in range(num_lays_height):
        conds = ['maxy', ]

        # expand directions (normal to border)
        expand_dirs = [
            np.array([0, 1]),
        ]

        # start-end direction (parallel to border)
        end_dirs = [
            np.array([1, 0]),
        ]

        nd, ndlab, el, ellab = expand_mult_dirs(
            conds, expand_dirs, end_dirs, nd, ndlab, el, ellab,
            num_layers=[1, ], width=width_lays_height
        )

    # ****** ADD LOADING HOLE *******
    circ_centre = (
        - (specimen['W'] - specimen['A'] / 2),
        specimen['F'] / 2
    )
    circ_rad = specimen['C'] / 2

    pts_on_circ = utils.circle_points(circ_rad, 80, circ_centre)
    pts_on_circ_ridge1 = utils.circle_points(circ_rad - 2, 80, circ_centre)

    # FIND NODES CLOSE TO HOLE
    atol = 13
    nd_dist_circ = (nd - circ_centre)[:, 0] ** 2 + (nd - circ_centre)[:, 1] ** 2
    pts_in_idx = np.nonzero(nd_dist_circ < (circ_rad - 0.2) ** 2)
    pts_out_idx = np.nonzero(nd_dist_circ >= (circ_rad - 0.2) ** 2)[0]
    pts_close_idx = np.nonzero(np.isclose(
        np.abs(nd_dist_circ - (circ_rad - 0.2) ** 2), 0, atol=atol))
    pts_out_close_idx = np.nonzero(
        np.logical_and(
            np.isclose(np.abs(nd_dist_circ - (circ_rad - 0.2) ** 2), 0, atol=atol),
            #         (pts_close_idx, nd_dist_circ >= (circ_rad - 0.2) ** 2 - 0.2)[1]
            nd_dist_circ >= (circ_rad - 0.2) ** 2 - 0.2
        )
    )[0]

    # nodes hole top
    ridge_el = []
    ridge_nd_top = []
    ridge_ndlab_top = []

    # loop through points near circle
    for i, pt in enumerate(nd[pts_out_close_idx]):
        # distance between pt and pts on circle
        diff_vecs = pt - pts_on_circ
        dist = (diff_vecs[:, 0] ** 2 + diff_vecs[:, 1] ** 2) ** 0.5
        nd[pts_out_close_idx[i]] = pts_on_circ[np.argmin(dist)]

        if pts_on_circ[np.argmin(dist)][1] > 22.0:
            ridge_nd_top.append(pts_on_circ[np.argmin(dist)][None])
            ridge_ndlab_top.append(ndlab[pts_out_close_idx[i]])

    ridge_sort_idx = np.argsort(np.concatenate(ridge_nd_top, axis=0)[:, 0])
    ridge_nd_top_sort = np.array(ridge_ndlab_top)[ridge_sort_idx]
    for i, n in enumerate(ridge_nd_top_sort[:-1]):
        ridge_el.append([
            np.max(ndlab)+1, np.max(ndlab)+1, ridge_nd_top_sort[i+1], n
        ])
    ridge_el = np.array(ridge_el)

    # REMOVE ELEMENTS
    el_h = np.copy(el)
    ellab_h = np.copy(ellab)

    for i, e in enumerate(el_h):
        if np.any(np.isin(e, ndlab[pts_in_idx])):
            #         print(np.nonzero(np.isin(e, ndlab[pts_in_idx]))[0])
            #         print(ellab_h[i], e)
            rem_idx = np.nonzero(np.isin(e, ndlab[pts_in_idx]))[0]
    #         print('rem_idx: ', rem_idx)
            if rem_idx.shape[0] == 1:
                if rem_idx[0] < len(e)-1:
                    e[rem_idx[0]] = e[rem_idx[0]+1]
                else:
                    e[rem_idx[0]] = e[0]
                el_h[i] = e

    rem_idx = np.concatenate((np.nonzero(np.count_nonzero(np.isin(el_h, ndlab[pts_in_idx]), axis=1) == 2)[0],
                              np.nonzero(np.count_nonzero(
                                  np.isin(el_h, ndlab[pts_in_idx]), axis=1) == 3)[0],
                              ))
    rem_idx = np.concatenate((rem_idx, np.nonzero(
        np.all((np.isin(el, ndlab[pts_in_idx])), axis=1))[0]))

    el_h = np.delete(el_h, rem_idx, axis=0)
    ellab_h = np.delete(ellab_h, rem_idx, axis=0)

    # CREATE RIDGE
    nd = nd[pts_out_idx]
    ndlab = ndlab[pts_out_idx]
    nd = np.append(nd, np.array(circ_centre)[None], axis=0)
    # print(nd)
    ndlab = np.append(ndlab, np.max(ndlab)+1)
    ridge_centre_lab = ndlab[-1]

    specimen_set = (np.min(ellab_h), np.max(ellab_h))
    el_h = np.concatenate((el_h, ridge_el))

    ridge_el_set = np.arange(np.max(ellab_h)+1, np.max(ellab_h)+1+ridge_el.shape[0])
    ellab_h = np.concatenate((ellab_h, ridge_el_set))

    num_ridge_el = ridge_el.shape[0]
    num_nodes_layer = len(ndlab)
    loadline_nset = ndlab[-1]
    mid_plane_set = ndlab

    # **** Find BCs ****
    nodes_b, labs_b = find_border(nd, ndlab,
                                  condition='miny')

    out = {
        'node_coordinates': nd,
        'node_labels': ndlab,
        'element_nodes': el_h,
        'element_labels': ellab_h,
        'elsets': {
            'specimen': specimen_set,  # defined as list or (start, stop, [step])
            'ridge': ridge_el_set,

        },
        'nsets': {
            'crackfront': labs_b,
            'load-line': loadline_nset,
            #         'cracktip': crack_nodes['crack_line'][:,0],
            #         'crackline': crack_nodes['crack_line'][:,1:].flatten(),
            'cracklip': crack_nodes['crack_lip'],
            'flank': np.concatenate((crack_nodes['crack_line'].flatten(), crack_nodes['crack_lip'])),
            'midplane': mid_plane_set,
        },
    }

    return out
