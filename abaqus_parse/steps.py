def generate_compact_tension_specimen_steps(applied_displacement, number_contours):
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
            'time increments': (0.02, 1.0, 1e-08, 0.02),
            'bcs':[{
                'node set': 'load-line',
                'dof': (2, 2, applied_displacement)
            },],
            'output':
            {
                'restart frequency': 0,
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
                'field':{},
            }
            
        }
    }

    return steps