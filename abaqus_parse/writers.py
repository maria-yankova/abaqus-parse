import numpy as np
from abaqus_parse.utils import format_arr 
__all__ = []

def write_inp(path, geometry, sections, materials, steps):
    """
    Write the input file '.inp' for an Abaqus simulation.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    """
    # ********** GEOMETRY **********
    nodes, node_labs, elems, elem_labs, elem_type, elem_sets, node_sets = geometry.values()
    
    sep = '\n'

    nodes = [
        '**\n',
        '**Nodes\n',
        '**\n',
        '*Node\n',
        format_arr([node_labs[None].T, nodes],
                            format_spec=['{:d}', '{:12.7f}'],
                            col_delim=', ')  
    ]
    
    elems = [
        '**\n',
        '**Elements\n',
        '**\n',
        '*Element, type=' + elem_type + '\n',
        format_arr([elem_labs[None].T, elems],
                            format_spec=['{:d}', '{:d}'],
                            col_delim=', ')
    ]
    
    n_sets = [
        '**\n',
        '**Node sets\n',
        '**\n',
    ]
    for k, v in node_sets.items():
        if k=='load-line':
            print(v, type(v))
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
        elif type(v)==np.int32:
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

        
            
    geom = sep.join([
        ''.join(nodes),
        ''.join(elems),
        ''.join(el_sets),
    ])
    
    # ********** SECTIONS & MATERIALS **********
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
    mats = [
        '**\n',
        '**Materials\n',
        '**\n',
    ]
    for k, v in materials.items():
        mats.append(
            '*Material, name=' + k + '\n')
        for sk, ss in v.items():
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
                '*' + v['type'] + '\n' + format_arr(list(np.array(list(v['time increments']))[None].T), format_spec=['{:3.2f}', '{:3.1f}','{:2.1e}', '{:3.2f}'],
                                                   col_delim=', ')
            )
        for bc in v['bcs']:
                # print(list(np.array((bc['dof']))[None].T))
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
                    if 'cracks' in vo.keys():
                        for crack in vo['cracks']:
                            stps.append(
                                '*Contour Integral, crack name=' + str(crack['name']) + ', contours=' + str(crack['contours'])
                            )
                            if 'crack tip nodes' in crack.keys():
                                stps.append(', crack tip nodes')
                                if crack['symmetry']:
                                    stps.append(', symm')
                                stps.append('\n')
                                if any(isinstance(el, list) for el in crack['crack tip nodes']):
                                    for cr in crack['crack tip nodes']:
                                         stps.append(cr[0] + ', ' + cr[1] + ', ' + format_arr(np.array(crack['direction']), format_spec=['{:d}'], col_delim=', '))
                                else:
                                    stps.append('\n' + crack['crack tip nodes'][0] + ', ' + crack['crack tip nodes'][1] + ', ' + format_arr(np.array(crack['direction']), format_spec=['{:d}'], col_delim=', '))
                            
        if k != 'initial-step':       
            stps.append('*End Step\n')

#     all_docs = ''.join([
#         ''.join(nodes),
#         ''.join(elems),
#         ''.join(el_sets),
#         ''.join(sects),
#         ''.join(mats),
#         ''.join(stps),
#     ])
    
    with open(path, 'w') as of:   # to do: test if i need universal newline mode here

        # Input file heading
        of.write('*Heading\n')
        of.write('** Job name: job-91 Model name: Model-1\n')
        of.write('** Generated by: Abaqus/CAE 2017\n')
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

   