import numpy as np
from abaqus_parse.utils import format_arr
from abaqus_parse._version import __version__
__all__ = []

def write_inp(path, materials, parts, steps, assembly=None):
    """
    Write the input file '.inp' for an Abaqus simulation.
    
    Parameters
    ----------
    path : string
        Path for writing the input file.
    materials : dict of dicts
         Dict whose keys are the material names and whose values are dicts with the following keys:
            elastic : dict, optional
                Dict with keys:
                    youngs_modulus : float
                    poisson_ratio : float
            func_name : dict, optional, a function that generates a dictionary with the material properties, e.g. 
                elastic :  dict
                    Dict with keys:
                        youngs_modulus : float
                        poisson_ratio : float
                plastic : dict
                    stress_strain : ndarray of shape (R, 2)
            **func_kwargs, optional, input arguments for func_name, e.g. temperature.
    parts : dict of dicts
        Dict whose keys are the part names and whose values are dicts with the following keys:
            node_coordinates : ndarray of shape (M, N) for M nodes with N coordinates
            node_labels : ndarray of shape (M,) for M node labels
            element_nodes : ndarray of shape (L, K) for L elements with K nodes each
            element_labels : ndarray of shape (L,) for L label nodes
            nsets : dict
                Dict keys are the node set names and values are one of:
                    - range of nodes specified as a tuple (start, stop, [step]), OR
                    - ndarray of shape (P,) for P nodes, OR
                    - integer for a single node
            elsets : dict
                Dict keys are the element set names and values are one of:
                    - range of elements specified as a tuple (start, stop, [step]), OR
                    - ndarray of shape (Q,) for Q elements, OR
                    - integer for a single element
            sections : list of dicts
                Each dictionary has the following keys and values:
                    type : str
                        The type of section, e.g. 'Solid'.
                    material : str
                        The name of the material of the section.
                    elset : str
                        The name of the element set part of the section. 
    assembly : list of dicts
        Definitions of the part instances that the assembly consists of (optional) (TODO).
    steps: dict of dicts (to be changed to list of dicts)
        Dict whose keys are the step names and whose values are dicts with the following keys:
        'initial-step' : dict, compulsory
            Dictionary with a key 'bcs' and value a list of dictionaries with the following keys and values:
                node set : str 
                    The node set name
                dof : tuple, optional 
                    Degrees of freedom to be included into the boundary condition
                type : str, optional
                    Type of boundary condition, e.g. 'XSYMM' OR 'PINNED'.
  
    Returns
    -------
    An Abaqus .inp file.

    TODO:
    - Add preprint options as input parameters. Currently, hard coded.
    - User specified Field Output. Currently, hard coded.

    """
    # ********** PARTS **********

    for part, part_definition in parts.items():

        node_coordinates = part_definition['node_coordinates']
        node_labels = part_definition['node_labels']
        elem_nodes = part_definition['element_nodes']
        elem_labs = part_definition['element_labels']
        elem_type = part_definition['element_type']
        elem_sets = part_definition['elsets']
        node_sets = part_definition['nsets']
        sections = part_definition['sections']

        sep = '\n'

        nodes = [
            '**\n',
            '**Nodes\n',
            '**\n',
            '*Node\n',
            format_arr([node_labels[None].T, node_coordinates],
                                format_spec=['{:d}', '{:12.7f}'],
                                col_delim=', ')  
        ]
        
        elems = [
            '**\n',
            '**Elements\n',
            '**\n',
            '*Element, type=' + elem_type + '\n',
            format_arr([elem_labs[None].T, elem_nodes],
                                format_spec=['{:d}', '{:d}'],
                                col_delim=', ')
        ]
        
        n_sets = [
            '**\n',
            '**Node sets\n',
            '**\n',
        ]
        for k, v in node_sets.items():
            if type(v) == tuple:
                n_sets.append(
                    '*Nset, nset=' + k + ', generate\n' +\
                    str(v[0]) + ', ' + str(v[1]) + ', 1\n' 
                            )
            elif type(v)==list or type(v)==np.ndarray:
                if type(v)==list:
                    v = np.array(v)
                whole_rows = v.size // 16
                first_block = v[:(whole_rows * 16)].reshape(-1, 16)
                remaining_block = v[(whole_rows * 16):]
                n_sets.append('*Nset, nset=' + k + '\n' +\
                            format_arr(first_block, format_spec=['{:d}'],
                                    col_delim=', ') +\
                            format_arr(remaining_block, format_spec=['{:d}'],
                                    col_delim=', '))
            elif type(v)==np.int32 or type(v)==np.int64:
                n_sets.append('*Nset, nset=' + k + '\n' + str(v) + '\n') 
        el_sets = [
            '**\n',
            '**Element sets\n',
            '**\n',
        ]
        for k, v in elem_sets.items():
            if type(v) == tuple:
                el_sets.append(
                    '*Elset, elset=' + k + ', generate\n' +\
                    str(v[0]) + ', ' + str(v[1]) + ', 1\n' 
                            )
            elif type(v)==list or type(v)==np.ndarray:
                whole_rows = v.size // 16
                first_block = v[:(whole_rows * 16)].reshape(-1, 16)
                remaining_block = v[(whole_rows * 16):]
                n_sets.append('*Elset, elset=' + k + '\n' +\
                            format_arr(first_block, format_spec=['{:d}'],
                                    col_delim=', ') +\
                            format_arr(remaining_block, format_spec=['{:d}'],
                                    col_delim=', '))
            elif type(v)==np.int32 or type(v)==np.int64:
                n_sets.append('*Elset, elset=' + k + '\n' + str(v) + '\n') 
        
        # Sections
        sects = [
            '**\n',
            '**Sections\n',
            '**\n',
        ]
        for sect in sections:   
            sects.append(
                '*' + sect['type'] + ' Section, elset=' + sect['elset']  +\
                ', material=' + sect['material'] + '\n'
            )
        sects = sep.join([
            ''.join(sects),
        ])
     # ********** MATERIALS **********
    mats = [
        '**\n',
        '**Materials\n',
        '**\n',
    ]
    for k, v in materials.items():
        mats.append(
            '*Material, name=' + k + '\n')
        for sk, ss in v.items():
            if sk == 'Elastic':
                ss = [[ss['youngs_modulus'], ss['poisson_ratio']]]
            if sk == 'Plastic':
                ss = ss['stress_strain']
            mats.append(
                '*' + sk + '\n' +\
                format_arr(np.array(ss), format_spec=['{:12.7f}'],
                            col_delim=', ')
            )
            
    # ********** STEPS **********       
    stps = [
        '**\n',
        '**Boundary conditions\n',
        '**\n',
    ]
    for k, v in steps.items():
        if k != 'initial-step':
            
            stps.append(
                '*Step, name=' + v['name'] + ', nlgeom=YES\n' +\
                '*' + v['type'] + '\n' + format_arr(list(np.array(list(v['time_increment_definition']))[None].T), format_spec=['{:3.2f}', '{:3.1f}','{:2.1e}', '{:3.2f}'],
                                                   col_delim=', ')
            )
        for bc in v['bcs']:
                stps.append(
                    '*Boundary\n' + bc['node set'] +', ')
                if 'dof' in bc.keys():
                    if len(bc['dof'])==2:
                        stps.append(format_arr(np.array(bc['dof']), format_spec=['{:d}'], col_delim=', ') )
                    else:
                        stps.append(str(bc['dof'][0]) + ', ' + str(bc['dof'][1])+ ', ' + str(bc['dof'][2]) + '\n') #list(np.array((bc['dof']))[None].T)
                        # (format_arr(list(np.array(i) for i in bc['dof']), format_spec=['{:d}','{:d}', '{:3.1f}'], col_delim=', ') 
                elif 'type' in bc.keys():
                    stps.append(bc['type'] + '\n')
        if 'output' in v.keys():
            stps.append(
                '*Restart, write, frequency=' + str(v['output']['restart frequency']) +'\n'
            )
            stps.append(
                        '*Output, field\n*Node Output\nCOORD, U\n*Element Output, position=centroidal, elset=specimen\nE, EVOL, PE, PEEQ, S, COORD\n'
                    )
            for ko, vo in v['output'].items():
                
                if ko=='history':
                    stps.append(
                        '*Output, history, frequency=' + str(vo['frequency']) + '\n'
                    )
                    # if 'cracks' in vo.keys():
                    #     for crack in vo['cracks']:
                    #         stps.append(
                    #             '*Contour Integral, crack name=' + str(crack['name']) + ', contours=' + str(crack['contours'])
                    #         )
                    #         if 'crack tip nodes' in crack.keys():
                    #             stps.append(', crack tip nodes')
                    #             if crack['symmetry']:
                    #                 stps.append(', symm')
                    #             stps.append('\n')
                    #             if any(isinstance(el, list) for el in crack['crack tip nodes']):
                    #                 for cr in crack['crack tip nodes']:
                    #                      stps.append(cr[0] + ', ' + cr[1] + ', ' + format_arr(np.array(crack['direction']), format_spec=['{:d}'], col_delim=', '))
                    #             else:
                    #                 stps.append('\n' + crack['crack tip nodes'][0] + ', ' + crack['crack tip nodes'][1] + ', ' + format_arr(np.array(crack['direction']), format_spec=['{:d}'], col_delim=', '))
                            
        if k != 'initial-step':       
            stps.append('*End Step\n')
    
    with open(path, 'w') as of:   # to do: test if i need universal newline mode here

        # Input file heading
        of.write('*Heading\n')
        of.write('** Generated by: abaqus-parse v' + str(__version__) +' \n')
        of.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
        of.write(
            ''.join([
                ''.join(nodes),
                ''.join(elems),
                ''.join(n_sets),
                ''.join(el_sets),
                ''.join(sects),
                ''.join(mats),
                ''.join(stps),
            ])
        )

   