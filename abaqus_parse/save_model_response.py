import numpy as np

def save_model_response(fname):
    
    data = np.loadtxt(fname, comments='%')

    model_response = {
        'time': data[:, 0],
        'S_mid_mises': data[:, 1],
        'LE_mid_mises': data[:, 2],
        'LE_mid_11': data[:, 3],
        'LE_mid_22': data[:, 4],
        'S_corner_mises': data[:, 5],
        'LE_corner_mises': data[:, 6],
        'LE_corner_11': data[:, 7],
        'LE_corner_22': data[:, 8],
    }

    return model_response
