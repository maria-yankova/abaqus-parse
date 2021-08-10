import numpy as np
from warnings import warn


def compute_forming_limit_curve(all_model_responses, strain_rate_ratio_threshold,
                                all_displacement_BCs, all_groove_angles):
    """
    Parameters
    ----------
    all_model_responses : list of dict
        Each list item must be a dict with keys:
            time
            LE_mid_mises
            LE_corner_mises
            LE_corner_22
            LE_corner_11
    strain_rate_ratio_threshold : number
    all_displacement_BCs : list of list
    all_groove_angles : list of float

    Returns
    -------
    forming_limit_curve : dict
        Dict with the following keys:
            strain_paths : list of list of ndarray
                Each list element is a sub-list corresponding to the strain paths of a 
                given groove angle. Each sub-list item is an array with two rows (minor
                and major strain respectively), representing the strain path up to
                necking. 
            forming_limits : ndarray
                Each column is the most conservative forming limit (minor and major
                strain) for a given strain path.

    """

    forming_limit_curve = {
        'displacement_BCs': all_displacement_BCs,
        'groove_angles': all_groove_angles,
    }

    strain_paths = []

    for resp, groove_angle, disp_BCs in zip(
        all_model_responses,
        all_groove_angles,
        all_displacement_BCs
    ):

        first_derivative = np.diff(resp['LE_11'], axis=0)

        # For all increments, identify elements with min/max strain rate:
        elem_idx_min = np.argmin(first_derivative, axis=1)
        elem_idx_max = np.argmax(first_derivative, axis=1)

        # For all increments, find min/max strain rate over all elements
        min_derivative = first_derivative[:, elem_idx_min]
        max_derivative = first_derivative[:, elem_idx_max]

        # For all increments, maximal strain rate ratio
        min_derivative[np.isclose(min_derivative, 0)] = np.nan
        strain_rate_ratio = np.abs(max_derivative / min_derivative)

        strain_path_i = {
            'groove_angle_deg': groove_angle,
            'displacement_BCs': disp_BCs,
            'strain_rate_ratio': strain_rate_ratio,
            'strain': None,
        }
        strain_paths.append(strain_path_i)

        try:
            first_threshold_idx = np.where(
                strain_rate_ratio >= strain_rate_ratio_threshold
            )[0][0]
        except IndexError:
            warn(f'Simulation did not reach the target threshold value '
                 f'({strain_rate_ratio_threshold}) for groove angle '
                 f'{groove_angle} (degrees) and displacment BCs: {disp_BCs}.')
            continue

        idx_min = elem_idx_min[first_threshold_idx]

        first_threshold_idx = first_threshold_idx + 1

        minor_strain = resp['LE_22'][:first_threshold_idx + 1, idx_min]
        major_strain = resp['LE_11'][:first_threshold_idx + 1, idx_min]
        strain_path = np.vstack((minor_strain, major_strain))
        strain_paths[-1].update({'strain': strain_path})

    forming_limit_curve.update({'strain_paths': strain_paths})

    # Get the overal forming limit for each strain path, taking the smallest major strain
    # over all groove angles for a given strain path:
    forming_limits = []
    seen_disp_BCs = []
    selected_angles = []
    for disp_BC in all_displacement_BCs:

        if disp_BC not in seen_disp_BCs:
            seen_disp_BCs.append(disp_BC)
        else:
            continue

        strain_paths_sub = [i for i in strain_paths if i['displacement_BCs'] == disp_BC]

        if all([i['strain'] is None for i in strain_paths_sub]):
            warn(f'No simulations (at any groove angle) reached the target threshold '
                 f'value for disp_BC: {disp_BC}.')
            forming_limits.append([np.nan, np.nan])
            selected_angles.append(np.nan)
            continue

        all_max_maj = []
        # Looping over groove angles:
        for strain_path_i in strain_paths_sub:
            if strain_path_i['strain'] is None:
                all_max_maj.append(np.inf)
            else:
                all_max_maj.append(np.max(strain_path_i['strain'][1]))

        safest_angle_idx = np.argmin(all_max_maj)
        forming_limits.append(strain_paths_sub[safest_angle_idx]['strain'][:, -1])
        selected_angles.append(strain_paths_sub[safest_angle_idx]['groove_angle_deg'])

    forming_limits = np.array(forming_limits)
    selected_angles = np.array(selected_angles)
    srt_idx = np.argsort(forming_limits[:, 0])  # Sort by minor strain
    final_forming_limits = forming_limits[srt_idx].T  # Transpose to column vectors
    selected_angles = selected_angles[srt_idx]

    forming_limit_curve.update({
        'forming_limits': final_forming_limits,
        'forming_limit_groove_angles_deg': selected_angles,
    })

    return forming_limit_curve
