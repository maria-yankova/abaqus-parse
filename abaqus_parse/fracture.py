def standard_specimen(spec_type, dimensions='2D', fraction='half', aw_ratio=0.5):
    """
    Returns the dimensions of standard fracture specimens.

    Parameters
    ----------
    spec_type
        CT-1T/SENB
    dimensions
        2D/3D
    fraction
        quarter/half/full
    aw_ratio
        fix at 0.5 for CT
    
    Returns
    -------
    if 2D:
        a/w, W, A, C, D, E, F
    if 3D:h
        a/w, W, A, B, C, D, E, F
    """
    
    specimens_dims = {
    'ct-1t':{
            'a/w': 0.5,
            'W': 50,
            'A': 62.5,
            'B': 25,
            'C': 12.5,
            'D': 23.5,
            'E': 60,
            'F': 37.5
        }
    }
    specimen = specimens_dims[spec_type]
    if dimensions=='2D':
        specimen.pop('B')
        
    if fraction=='half':
        specimen['E'] *= 0.5

    return specimens_dims['ct-1t']
    