import numpy as np


def save_model_response(fname):

    data = np.loadtxt(fname, comments='%')

    nb_el = int(np.where(data[:, 0] > 0.0)[0][0])
    nb_frame = int(np.shape(data)[0]/nb_el)

    time = data[0:-1:nb_el, 0]
    LE_mises = np.reshape(data[:, 2], [nb_frame, nb_el])
    LE_11 = np.reshape(data[:, 3], [nb_frame, nb_el])
    LE_22 = np.reshape(data[:, 4], [nb_frame, nb_el])

    model_response = {
        'time': time,
        'LE_mises': LE_mises,
        'LE_11': LE_11,
        'LE_22': LE_22,
    }

    return model_response
