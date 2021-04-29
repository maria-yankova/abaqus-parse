import numpy as np
from warnings import warn


def compute_forming_limit_curve(all_model_responses, strain_rate_ratio_threshold,
                                num_groove_angles):
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
    num_groove_angles : int

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

    num_strain_paths = int(len(all_model_responses) / num_groove_angles)
    forming_limit_curve = {
        'strain_paths': [[None for _ in range(num_groove_angles)]
                         for _ in range(num_strain_paths)],
        'forming_limits': None,
    }

    for resp_idx, resp in enumerate(all_model_responses):

        groove_angle_idx = resp_idx // num_strain_paths
        strain_path_idx = resp_idx % num_strain_paths

        dt = np.gradient(resp['time'])
        first_derivative_mid = np.gradient(resp['LE_mid_mises']) / dt
        first_derivative_corner = np.gradient(resp['LE_corner_mises']) / dt

        strain_rate_ratio = first_derivative_mid / first_derivative_corner

        # loop over all corners?
        try:
            first_threshold_idx = np.where(
                strain_rate_ratio >= strain_rate_ratio_threshold
            )[0][0]
        except IndexError:
            warn(f'Simulation did not reach the target threshold value '
                 f'({strain_rate_ratio_threshold}) for groove_angle_idx '
                 f'{groove_angle_idx} and strain_path_idx: {strain_path_idx}.')
            continue

        minor_strain = resp['LE_corner_22'][:first_threshold_idx + 1]
        major_strain = resp['LE_corner_11'][:first_threshold_idx + 1]
        strain_path = np.vstack((minor_strain, major_strain))
        forming_limit_curve['strain_paths'][strain_path_idx][groove_angle_idx] = strain_path

    # Get the overal forming limit for each strain path, taking the smallest major strain
    # over all groove angles for a given strain path:
    forming_limits = []
    for strain_path_idx in range(num_strain_paths):

        if all([i is None for i in forming_limit_curve['strain_paths'][strain_path_idx]]):
            warn(f'No simulations (at any groove angle) reached the target threshold '
                 f'value for strain_path_idx: {strain_path_idx}.')
            forming_limits.append([np.nan, np.nan])
            continue

        all_max_maj = []
        for groove_angle_i_path in forming_limit_curve['strain_paths'][strain_path_idx]:
            if groove_angle_i_path is None:
                all_max_maj.append(np.inf)
            else:
                all_max_maj.append(np.max(groove_angle_i_path[1]))

        safest_angle_idx = np.argmin(all_max_maj)
        forming_limits.append(
            forming_limit_curve['strain_paths'][strain_path_idx][safest_angle_idx][:, -1]
        )

    forming_limit_curve['forming_limits'] = np.array(forming_limits).T

    return forming_limit_curve
