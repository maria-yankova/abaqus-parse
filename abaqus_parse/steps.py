def generate_compact_tension_specimen_steps(applied_displacement, number_contours, time_increment_definition):

    time_increment_def = (
        time_increment_definition['initial_time_increment'],
        time_increment_definition['total_step_time'],
        time_increment_definition['min_time_increment_allowed'],
        time_increment_definition['max_time_increment_allowed'],
    )
    steps = {
        'initial-step': {
            'bcs': [
                {
                'node set': 'load-line',
                'dof': (1, 1)
                },
                {
                'node set': 'load-line',
                'dof': (2, 2)
                },
                {
                    'node set': 'crackfront',
                    'type': 'YSYMM'
                },
                {
                    'node set': 'midplane',
                    'type': 'ZSYMM'
                }
            ]
        },
        'load-step-1': {
            'name': 'Step-1',
            'type': 'Static',
            'time_increment_definition': time_increment_def, #(0.02, 1.0, 1e-08, 0.02),
            'bcs':[{
                'node set': 'load-line',
                'dof': (2, 2, applied_displacement)
            },],
            'output':
            {
                'restart frequency': 0,
                'field':[
                    {
                        'output type': 'node',
                        'variables': ['COORD', 'U'],
                    },
                    {
                        'output type': 'element',
                        'position': 'centroidal',
                        'elset': 'specimen',
                        'frequency': 0,
                        'variables': ['E', 'EVOL', 'PE', 'PEEQ', 'S', 'COORD'],
                    }
                ],
                'history':
                {
                    'frequency': 0,
                    'cracks': [
                        {
                            'name': 'CRACK',
                            'contours': number_contours,
    #                         'crack tip nodes': [ 'crackline', 'cracktip'],
                            'crack tip nodes': [],
                            'symmetry': True,
                            'direction': [1, 0, 0],
                        }
                    ]
                },
            }
            
        }
    }

    return steps