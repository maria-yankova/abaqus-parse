import numpy as np

def get_euromat_A_material_props(temp):
    """Get mechanical properties of Euromaterial A, calculated as a function of
    temperature from a fit to uniaxial tests [1].
    
    Parameters
    ----------
    temp : float
        Temperature in Celsius.
    
    References
    ---------
    [1] Beardsmore, D. & Lidbury, D. (2006). Serco Assurance report. PERFECT.
    
    Returns
    -------
    material_props : dict
        Dict containing keys:
            yield_stress : float
                Yield stress in MPa.
            youngs_modulus : float
                Young's modulus in MPa.
            hardening_exponent : float
            UTS : float
                Ultimate tensile strength in MPa.
            poisson_ratio : float  
    """
    
    if temp >= 0:
        hard_exp = 9.96888 - 0.00185 * temp
    else:
        hard_exp = (9.96888 - 901.98738 * temp - 3.01593 * temp ** 2) / (1 - 90.48011 * temp - 0.41856 * temp ** 2)
    
    material_props = {
        'yield_stress': 421.18 + 63.9 * np.exp(-temp / 91),
        'youngs_modulus': -90 * temp + 206000,
        'hardening_exponent': hard_exp,
        'UTS': 564.1 + 70.2 * np.exp(-temp / 108),
        'poisson_ratio': 0.3,
    }
    
    return material_props


def get_euromat_A_material_definition(temp, plastic_stress_max=1200, plastic_num_incs=200):
    """
    Get an Abaqus compatible material definition for Euromaterial A.
    
    Paramaters
    ----------
    temp : float
        Temperature in Celsius.
    plastic_stress_max : float or int
        Maximum plastic stress in MPa.
    plastic_num_incs : int
        Number of increments of the plastic stress-strain return data.
    
    Returns
    -------
    material_defn : dict
        Dict containing keys:
            elastic : dict
                Dict containing keys:
                    youngs_modulus : float
                        Young's modulus in MPa.
                    poisson_ratio : float
            plastic : dict
                Dict containing keys:
                    stress_strain : ndarray of shape (N, 2)
                        Plastic stress-strain data; first column is stress.
    """
    
    euromat_props = get_euromat_A_material_props(temp)

    s0_wt = euromat_props['yield_stress']
    E_wt = euromat_props['youngs_modulus']
    n_wt = euromat_props['hardening_exponent']
    
    stress_data = np.linspace(s0_wt, plastic_stress_max, plastic_num_incs)
    tot_train = (s0_wt / E_wt) * (stress_data / s0_wt) ** n_wt
    pl_strain_data = tot_train - stress_data / E_wt
    
    material_defn =  {
        'Elastic': {
            'youngs_modulus' : E_wt,
            'poisson_ratio' : euromat_props['poisson_ratio'],
        }, 
        'Plastic': {
            'stress_strain': np.array([stress_data, pl_strain_data]).T
        }
    }
    return material_defn