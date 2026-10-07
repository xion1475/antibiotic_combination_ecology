# import packages
import cobra
import cometspy as c
import os
import pandas as pd
import re
import math
import itertools

# load models
E_WT_iM = cobra.io.read_sbml_model("../../models/iML1515.xml")
S0 = cobra.io.read_sbml_model("../../models/STM_v1_0_S0.xml")

# create E0 knockouts
E0_iM = E_WT_iM.copy()
E0_iM.genes.b3939.knock_out()

met_limiting_medium = {
    'EX_lcts_e': 10,
    'EX_met__L_e': .2,
    'EX_pi_e': 10,
    'EX_fe3_e': 10,
    'EX_mn2_e': 10,
    'EX_fe2_e': 10,
    'EX_zn2_e': 10,
    'EX_mg2_e': 10,
    'EX_ca2_e': 10,
    'EX_ni2_e': 10,
    'EX_cu2_e': 10,
    'EX_cobalt2_e': 10,
    'EX_mobd_e': 10,
    'EX_so4_e': 10,
    'EX_k_e': 10,
    'EX_cl_e': 10,
    'EX_o2_e': 10,
    'EX_nh4_e': 10}

carbon_limiting_medium = {
    'EX_gal_e': 3,
    'EX_met__L_e': 10,
    'EX_pi_e': 10,
    'EX_fe3_e': 10,
    'EX_mn2_e': 10,
    'EX_fe2_e': 10,
    'EX_zn2_e': 10,
    'EX_mg2_e': 10,
    'EX_ca2_e': 10,
    'EX_ni2_e': 10,
    'EX_cu2_e': 10,
    'EX_cobalt2_e': 10,
    'EX_mobd_e': 10,
    'EX_so4_e': 10,
    'EX_k_e': 10,
    'EX_cl_e': 10,
    'EX_o2_e': 10,
    'EX_nh4_e': 10}

# load functions
def get_gene_ids(model, genes):
    all_gene_names = [gene.name for gene in model.genes]
    all_gene_ids = [gene.id for gene in model.genes]
    names_and_ids = dict(zip(all_gene_ids, all_gene_names))
    gene_ids = [key for key, value in names_and_ids.items() if value in genes]
    return(gene_ids)

def scale_reaction(model: "Model", reaction_id: str, alpha: float, direction='forward'):
    if direction == 'forward':
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v >= 0)    # obtain 'end product'
    else:  # reverse reaction
        model.reactions.get_by_id(reaction_id).upper_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v <= 0)
    if isinstance(alpha, pd.Series):
        alpha = float(alpha)
    for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
        model.reactions.get_by_id(reaction_id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite

def separate_reaction(model: "Model", reaction_id: str, alpha: float):
    (lb, ub) = model.reactions.get_by_id(reaction_id).bounds
    rxn_ids = [reaction_id] #? 
    if(lb < 0 and ub !=0): # only perform if reaction is bidirectional
        rev_reaction = model.reactions.get_by_id(reaction_id).copy() # copy of target reaction
        rev_reaction.id = f'{reaction_id}_v1' # redefine id for copy of reaction
        model.add_reactions([rev_reaction]) # add to model
    
        scale_reaction(model, rev_reaction.id, alpha, direction='backward')
        
        rxn_ids.append(rev_reaction.id) # list of id for reaction and reaction_v1
    scale_reaction(model, reaction_id, alpha, direction='forward')
    return(rxn_ids)

# function to scale the stochiometry for all reactions related to a gene
def scale_stochiometry(model: "Model", alpha: float, gene: str, KO = False):
    gene_ids = get_gene_ids(model, gene)
    reactions = [reaction.id for reaction in model.genes.get_by_id(gene_ids[0]).reactions]
    rxn_ids = list()
    for rxn in reactions:
        if KO:
            model.reactions.get_by_id(rxn).knock_out()
        elif (rxn not in rxn_ids):
                rxn_ids.extend(separate_reaction(model, rxn, alpha))# copy of reaction, forward_terminal_change = True
    return(rxn_ids)

def get_alpha_value(model: "Model", start_alpha: float, gene: str, tolerance: float, frac_of_biomass: float, max_iterations: int = 100):
    test_model = model.copy()
    initial_biomass = test_model.slim_optimize() 
    try:
        value = float(initial_biomass)
    except ValueError:
        print("Non-float value detected. No initial biomass.")
        
    alpha = start_alpha
    low = 2
    high = None
    iteration = 0

    gene_id = get_gene_ids(test_model, gene)
    rxn_to_ko = [reaction.id for reaction in test_model.genes.get_by_id(gene_id[0]).reactions]
    knockout_model = test_model.copy()
    for rxn in rxn_to_ko:
        knockout_model.reactions.get_by_id(rxn).knock_out()
    essential = round(knockout_model.slim_optimize(), 5)

    if not math.isclose(essential, 0.0, abs_tol=1e-8):
        alpha = 100000
        
    else:
        while iteration < max_iterations:
            scale_stochiometry(test_model, alpha, gene)
            scaled_biomass = test_model.slim_optimize()
            
            if math.isclose(scaled_biomass, (initial_biomass * frac_of_biomass), rel_tol=tolerance):
                break
            
            if scaled_biomass > (initial_biomass * frac_of_biomass):  
                low = alpha
                alpha = alpha * 2
            else:
                high = alpha
                alpha = (low + high) / 2
            
            test_model = model.copy()
            iteration += 1
    
            if iteration == max_iterations:
                print(f"Warning: Maximum iterations ({max_iterations}) reached. Returning last tested alpha.")
                break

    return(alpha)

def make_scaled_model(model: "Model", alpha: float, gene: str):
    new_model = model.copy()
    scale_stochiometry(new_model, alpha, gene)
    return(new_model)

def pivot_long_fluxes(dataframe):
    df_long = pd.melt(
        dataframe,
        id_vars=['cycle', 'x', 'y'],
        var_name='reaction',
        value_name='flux'
    )
    
    if df_long[['x', 'y']].nunique().eq(1).all() and \
       (df_long['x'].iloc[0] == 1 and df_long['y'].iloc[0] == 1):
        df_long = df_long.drop(columns=['x', 'y'])

    df_long = df_long.drop_duplicates()

    return(df_long)

def fix_biomass_data(dataframe):
    
    new_name = 'biomass'  # or whatever you want to rename it to
    other_col = [col for col in dataframe.columns if col != 'cycle'][0]
    df = dataframe.rename(columns={other_col: new_name})

    return(df)

def drop_x_y(dataframe):
    if (dataframe['x'].eq(1).all()) and (dataframe['y'].eq(1).all()):
        df = dataframe.drop(columns=['x', 'y'])
    df = df.drop_duplicates()

    return(df)

def run_standard_monoculture_simulation(model: "Model", carbon: str, carbon_mmol: float, auxotrophies: str, auxo_mmol: float, max_cycles: int):
    comets_model = c.model(model)
    comets_model.open_exchanges()
    
    comets_model.initial_pop = [0, 0, 1.e-8]
    comets_model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
    
    test_tube = c.layout()
    test_tube.add_model(comets_model)
    
    base_nutrients = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e', 'h_e', 'k_e', 'h2o_e', 'mg2_e',
             'mn2_e', 'mobd_e', 'na1_e', 'ni2_e', 'nh4_e', 'o2_e', 'pi_e', 'so4_e', 'zn2_e']
    for nutrient in base_nutrients:
        test_tube.set_specific_metabolite(nutrient, 1000)
        test_tube.set_specific_static(nutrient, 1000)
    test_tube.set_specific_metabolite(carbon, carbon_mmol)
    if auxotrophies != "none":
        test_tube.set_specific_metabolite(auxotrophies, auxo_mmol)
    
    sim_params = c.params()
    sim_params.set_param("defaultKm", 0.00001) # M 
    sim_params.set_param("defaultVmax", 10) #mmol/gDw/hr
    sim_params.set_param("maxCycles", max_cycles)
    sim_params.set_param("timeStep", 1)
    sim_params.set_param("writeMediaLog", False)
    sim_params.set_param("writeFluxLog", False)
    
    comp_assay = c.comets(test_tube, sim_params)
    comp_assay.run()

    return(fix_biomass_data(comp_assay.total_biomass))

def run_standard_coculture_simulation(model1: "Model", model2: "Model", carbon: str, carbon_mmol: float, max_cycles: int):
    comets1 = c.model(model1)
    comets2 = c.model(model2)
    comets1.open_exchanges()
    comets2.open_exchanges()
    
    comets1.initial_pop = [0, 0, 1.e-8]
    comets2.initial_pop = [0, 0, 1.e-8]
    comets1.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" #this is for pFBA - MAX_OBJECTIVE_FLUX for fba
    comets2.obj_style =  "MAX_OBJECTIVE_MIN_TOTAL"
    
    test_tube = c.layout([comets1, comets2])
    
    base_nutrients = ["ca2_e", "cl_e", "cobalt2_e", "cu2_e","fe2_e", "fe3_e", "k_e","mg2_e",
          "mn2_e", "mobd_e", "ni2_e", "o2_e", "pi_e", "nh4_e", "so4_e", "zn2_e"]
    for nutrient in base_nutrients:
        test_tube.set_specific_metabolite(nutrient, 1000)
    test_tube.set_specific_metabolite(carbon, carbon_mmol)
    
    p = c.params()
    p.set_param("defaultKm", 0.00001) # M 
    p.set_param("defaultVmax", 10) #mmol/gDw/hr
    p.set_param("maxCycles", max_cycles)
    p.set_param("timeStep", 1)
    p.set_param("writeMediaLog", False)
    p.set_param("writeFluxLog", False)
    
    comp_assay = c.comets(test_tube, p)
    comp_assay.run()

    return(comp_assay.total_biomass)

def rename_ac_in_bulk(model, id):
    """Change name of acetate metabolite to a single bulk acetate
    
    Parameters
    ----------
    model: cobra.Model
    id: str
        id of metabolite to be renamed
    
    Returns
    -------
    None
    """
    
    new_id = f'bulk_{id}'
    model.metabolites.get_by_id(id).id = new_id
    if 'ac' in id:
        model.metabolites.get_by_id(new_id).name = 'bulk Acetate'
        model.metabolites.get_by_id(new_id).formula = 'C20H30O20'
    if '_e' in id:
        model.reactions.get_by_id(f'EX_{id}').id = f'EX_{new_id}'

def weight_carbon_byproduct(model, ac_scale, gal_scale):
    """Weight carbon byproduct secretion flux for E.coli with specific scales
    If a scale is None, no change is made
    Not changing functionality of metabolic model unlike stoichiometric scaling
    
    Parameters
    ----------
    model: Cobra.model
    all_components: dict
    ac_scale: float
        scale need to be larger than 1, reduce in flux accounted in pfba optimization procedure
    gal_scale: float
        scale need to be less than 1, increase in flux accounted in pfba optimization procedure

    Returns
    -------
    None
    """
    if gal_scale is not None:
        gal_scale = -1*(1-1/gal_scale) if gal_scale > 1 else -1*(1-gal_scale) # ensure increase in cose
    
    for rxn in [model.reactions.get_by_id('LACZpp')]:
        # TODO: only scale LACZpp is desirable, make gal secretion similar to knockout
        metab_to_scale = rxn.metabolites
        rxn.add_metabolites({k:v*gal_scale for k,v in metab_to_scale.items()}) # similar toE0.reactions.GALtex.knock_out()
    
    if ac_scale is not None: # flux do not count in minflux
        add_scale = ac_scale-1 if ac_scale>=1 else ac_scale
    
    # ac_p + 10 h_p <=> 10 ac_c + 10 h_c
    rxn = model.reactions.get_by_id("ACt2rpp")
    if not 'EX_' in rxn.id:
        metab_to_scale = rxn.metabolites
        rxn.add_metabolites({k:v*add_scale for k,v in metab_to_scale.items()})


E0 = E0_iM.copy()

ac_scale: float = 10
gal_scale: float = 3

add_metab_quantity_E = 1-1/ac_scale
add_metab_quantity_S = 1-1/ac_scale
E0.reactions.ACt2rpp.add_metabolites({E0.metabolites.ac_p: add_metab_quantity_E}) # 0.1 ac_p <-> ac_c
S0.reactions.ACt2rpp.add_metabolites({S0.metabolites.ac_p: add_metab_quantity_S})

weight_carbon_byproduct(E0, ac_scale, gal_scale)

for model, metabolite in itertools.product([E0, S0], ['ac_p', 'ac_e']):
    rename_ac_in_bulk(model, metabolite)

# single drug example
# load genes to knockdown
alpha_table = pd.read_csv("../../model-data/alpha-values/alpha-table.csv")
gene_list = alpha_table["gene"].to_list()

# single gene KOs
# single gene KOs for each strain
num_single_repression_cycles = 500
carbon_conc = 0.000278
met_conc = 0.000278 / 50

single_gene_biomasses = pd.DataFrame()

#base file
base = "../../model-data/alpha-values/single-drug/"

for i in range(0, len(gene_list)):
    gene = gene_list[i]

    # run monoculture S0 simulations
    S0_gene = alpha_table["S"][i]
    repressed_S0 = make_scaled_model(S0, S0_gene, gene)
    S_gene = run_standard_monoculture_simulation(model = repressed_S0, carbon = "gal_e", carbon_mmol = carbon_conc, auxotrophies = "none", auxo_mmol = 0, max_cycles = num_single_repression_cycles)

    # clean S monoculture dataframes
    S_gene["gene"] = gene
    S_gene["strain"] = "S"
    S_gene["interaction"] = "monoculture"
    
    # run monoculture E0 simulations
    E0_gene = alpha_table["E"][i]
    repressed_E0 = make_scaled_model(E0, E0_gene, gene)
    E_gene = run_standard_monoculture_simulation(model = repressed_E0, carbon = "lcts_e", carbon_mmol = carbon_conc, auxotrophies = "met__L_e", auxo_mmol = met_conc, max_cycles = num_single_repression_cycles)

    # clean S monoculture dataframes
    E_gene["gene"] = gene
    E_gene["strain"] = "E"
    E_gene["interaction"] = "monoculture"
        
    # run coculture simulation
    ES_gene = run_standard_coculture_simulation(model1 = repressed_E0, model2 = repressed_S0, carbon = "lcts_e", carbon_mmol = carbon_conc, max_cycles = num_single_repression_cycles)

    ES_gene = ES_gene.melt(
        id_vars=['cycle'],
        value_vars=[col for col in ES_gene.columns if col.startswith('iM') or col.startswith('STM')],
        var_name='strain',
        value_name='biomass'
    )

    ES_gene["gene"] = gene
    ES_gene["interaction"] = "coculture"

    # save individual data
    biomass = pd.concat([S_gene, E_gene, ES_gene])

    # concatenate data
    single_gene_biomasses = pd.concat([single_gene_biomasses, biomass])

    # print
    print("loop", i + 1, "of", len(gene_list) + 1)
    print("gene KO: ", gene)

single_gene_biomasses.to_csv(base + 'biomasses_unique_mono_no_gal.csv')