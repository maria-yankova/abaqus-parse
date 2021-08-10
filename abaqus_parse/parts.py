import numpy as np
from fracture_fea_laf.mesh import (
    make_donut_mesh,
    make_fine_plus_donut_mesh,
    compact_tension_specimen_mesh,
    bend_bar_specimen_mesh
)


def generate_compact_tension_specimen_parts(dimension, mesh_definition,
                                            elem_type, size_type, fraction, specimen_material):
    """
    Parameters
    -----------
    dimension: string
        Specification of the number of dimensions: '2D' or '3D'.
    """

    # refined_mesh_definition = make_fine_plus_donut_mesh(
    #     mesh_definition['crack_tip_radius_microns'],
    #     mesh_definition['fine_mesh_length'],
    #     mesh_definition['fine_mesh_element_length'],
    #     mesh_definition['fan_box_num_side_elements'],
    #     mesh_definition['fan_box_width'],
    #     ret_crack_definition=True,
    #     size_behind_crack=0.2
    # )
    # if dimension == '3D':
    #     refined_mesh_definition['number_layers'] = mesh_definition['number_layers']
    #     refined_mesh_definition['element_thickness'] = mesh_definition['element_thickness']

    return {
        'ct-specimen': {
            **compact_tension_specimen_mesh(mesh_definition, dimension, size_type, fraction),
            'element_type': elem_type,
            'sections':[
                {
                    'type': 'Solid',
                    'material': specimen_material,
                    'elset': 'specimen'    
                },
                {
                    'type': 'Solid',
                    'material': 'rigid',
                    'elset': 'ridge'
                }
            ],
        }
    }
    
def generate_bend_bar_specimen_parts(dimension, mesh_definition,
                     elem_type, size_type, fraction, specimen_material):

    """
    Parameters
    -----------
    dimension: string
        Specification of the number of dimensions: '2D' or '3D'.
    """

    refined_mesh_definition = make_fine_plus_donut_mesh(
        mesh_definition['crack_tip_radius_microns'],
        mesh_definition['fine_mesh_length'],
        mesh_definition['fine_mesh_element_length'],
        mesh_definition['fan_box_num_side_elements'],
        mesh_definition['fan_box_width'],
        ret_crack_definition=True,
        size_behind_crack=0.2
    )
    if dimension == '3D':
        refined_mesh_definition['number_layers'] = mesh_definition['number_layers']
        refined_mesh_definition['element_thickness'] = mesh_definition['element_thickness']

    return {
        'ct-specimen': {
            **bend_bar_specimen_mesh(refined_mesh_definition, dimension, size_type, fraction),
            'element_type': elem_type,
            'sections': [
                {
                    'type': 'Solid',
                    'material': specimen_material,
                    'elset': 'specimen'
                },
            ],
        }
    }
