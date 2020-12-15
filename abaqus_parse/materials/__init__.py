from abaqus_parse.materials import euromat_A

MATERIAL_MODEL_FUNCS = {
    'euromat_A': euromat_A.get_euromat_A_material_definition,
}

def generate_material_models(materials_list):
    
    material_models = {}
    
    for material in materials_list:
        mat_mod_i = {}

        if 'func_name' in material:
            func = MATERIAL_MODEL_FUNCS[material['func_name']]
            func_kwargs = {key: value for key, value in material.items()
                           if key not in ['name', 'func_name']}
            mat_mod_i = func(**func_kwargs) 

        else:
            if 'elastic' in material:
                mat_mod_i['Elastic'] = material['elastic']

        material_models[material['name']] = mat_mod_i


    return material_models