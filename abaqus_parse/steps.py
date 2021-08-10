def generate_compact_tension_specimen_steps(applied_displacement, number_contours, time_increment_definition, number_layers=None):

    time_increment_def = (
        time_increment_definition['initial_time_increment'],
        time_increment_definition['total_step_time'],
        time_increment_definition['min_time_increment_allowed'],
        time_increment_definition['max_time_increment_allowed'],
    )
    out = {
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
                    'node set': 'yplane',
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
            'time_increment_definition': time_increment_def,  # (0.02, 1.0, 1e-08, 0.02),
            'bcs': [{
                'node set': 'load-line',
                'dof': (2, 2, applied_displacement)
            }, ],
            'output':
            {
                'restart frequency': 0,
                'field': [
                    {
                        'output type': 'node',
                        'variables': ['COORD', 'U'],
                    },
                    {
                        'output type': 'element',
                        'position': 'centroidal',
                        'set name': 'specimen',
                        # 'frequency': 0,
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
                            # 'crack tip nodes': [ 'crackline0', 'cracktip0'],
                            'crack tip nodes': [],
                            'symmetry': True,
                            'direction': [1, 0, 0],
                        }
                    ]
                },
            }

        }
    }

    if number_layers:
        print('yes number_layers')
        for i in range(number_layers):
            out['load-step-1']['output']['history']['cracks'][0]['crack tip nodes'].append(
                ['crackline'+str(i), 'cracktip'+str(i)]
            )
    else:
        out['load-step-1']['output']['history']['cracks'][0]['crack tip nodes'].append(
            ['crackline', 'cracktip']
        )

    return out


def generate_threepoint_bend_test_steps(applied_displacement, number_contours, time_increment_definition, number_layers=None):

    time_increment_def = (
        time_increment_definition['initial_time_increment'],
        time_increment_definition['total_step_time'],
        time_increment_definition['min_time_increment_allowed'],
        time_increment_definition['max_time_increment_allowed'],
    )
    out = {
        'initial-step': {
            # 'interaction properties':[
            #      {
            #     'type': 'Surface Interaction',
            #     'dof': (1, 1)
            #     },
            # ]
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
                    'node set': 'load-line',
                    'dof': (3, 3)
                },
                {
                    'node set': 'load-line',
                    'dof': (4, 4)
                },
                {
                    'node set': 'load-line',
                    'dof': (5, 5)
                },
                {
                    'node set': 'load-line',
                    'dof': (6, 6)
                },
                {
                    'node set': 'yplane',
                    'type': 'YSYMM'
                },
                {
                    'node set': 'midplane',
                    'type': 'ZSYMM'
                },
                {
                    'node set': 'fixpin',
                    'type': 'ENCASTRE'
                }
            ]
        },
        'load-step-1': {
            'name': 'Step-1',
            'type': 'Static',
            'time_increment_definition': time_increment_def,  # (0.02, 1.0, 1e-08, 0.02),
            'bcs': [{
                'node set': 'load-line',
                'dof': (1, 1, applied_displacement)
            }, ],
            'output':
            {
                'restart frequency': 0,
                'field': [
                    {
                        'output type': 'node',
                        'variables': ['COORD', 'U'],
                    },
                    {
                        'output type': 'element',
                        'position': 'centroidal',
                        'set name': 'specimen',
                        # 'frequency': 0,
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
                            # 'crack tip nodes': [ 'crackline0', 'cracktip0'],
                            'crack tip nodes': [],
                            'symmetry': True,
                            'direction': [1, 0, 0],
                        }
                    ]
                },
            }

        }
    }

    if number_layers:

        for i in range(number_layers):
            out['load-step-1']['output']['history']['cracks'][0]['crack tip nodes'].append(
                ['crackline'+str(i), 'cracktip'+str(i)]
            )
    else:
        out['load-step-1']['output']['history']['cracks'][0]['crack tip nodes'].append(
            ['crackline', 'cracktip']
        )

    return out
