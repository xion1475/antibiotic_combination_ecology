# import packages
import cobra
import cometspy as c
import os
import pandas as pd
import re
import math

# load models
E_WT_iM = cobra.io.read_sbml_model("../../../models/iML1515.xml")
S0 = cobra.io.read_sbml_model("../../../models/STM_v1_0_S0.xml")

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
        print(f"Non-essential: returning alpha of 100000.")
    else:
        while iteration < max_iterations:
            scale_stochiometry(test_model, alpha, gene)
            scaled_biomass = test_model.slim_optimize()

            if math.isclose(scaled_biomass, initial_biomass * frac_of_biomass, rel_tol=tolerance):
                break

            if scaled_biomass > initial_biomass * frac_of_biomass:
                # Scaled too high — switch to binary search
                low = alpha
                if high is None:
                    alpha *= 2
                else:
                    alpha = (low + high) / 2
            else:
                # Scaled too low — record upper bound
                high = alpha
                alpha = (low + high) / 2

            test_model = model.copy()
            iteration += 1
            print(iteration)
            print(alpha)
    
            if iteration == max_iterations:
                print(f"Warning: Maximum iterations ({max_iterations}) reached. Returning last tested alpha.")
                break

            if round(alpha, 3) == 1:
                print(f"Warning: Any alpha > 1 reduces biomass > than the target threshold. Returning last alpha value.")
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
    ident = model.id
    
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
    sim_params.set_param("writeFluxLog", True)
    sim_params.set_param("FluxLogRate", 1)
    
    comp_assay = c.comets(test_tube, sim_params)
    comp_assay.run()

    return([fix_biomass_data(comp_assay.total_biomass), pivot_long_fluxes(comp_assay.fluxes_by_species[ident])])

def run_standard_coculture_simulation(model1: "Model", model2: "Model", carbon: str, carbon_mmol: float, max_cycles: int):
    id1 = model1.id
    id2 = model2.id
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
    p.set_param("writeFluxLog", True)
    p.set_param("FluxLogRate", 1)
    
    comp_assay = c.comets(test_tube, p)
    comp_assay.run()

    return([comp_assay.total_biomass, comp_assay.fluxes_by_species[id1], comp_assay.fluxes_by_species[id2]])

# single drug example
# load genes to knockdown
alpha_table = pd.read_csv("../../../model-data/alpha-values/alphas/alpha-table.csv")
gene_list = alpha_table["gene"].to_list()

# single gene KOs
# single gene KOs for each strain
num_single_repression_cycles = 750
carbon_conc = 0.000278
met_conc = 0.000278 / 50

single_gene_biomasses = pd.DataFrame()
single_gene_fluxes = pd.DataFrame()

#base file
base = "../../../model-data/alpha-values/single-drug/"

for i in range(0, len(gene_list)):
    gene = gene_list[i]

    # run monoculture S0 simulations
    S0_gene = alpha_table["S"][i]
    repressed_S0 = make_scaled_model(S0, S0_gene, gene)
    S_gene = run_standard_monoculture_simulation(model = repressed_S0, carbon = "gal_e", carbon_mmol = carbon_conc, auxotrophies = "none", auxo_mmol = 0, max_cycles = num_single_repression_cycles)

    # clean S monoculture dataframes
    S_gene[0]["gene"] = gene
    S_gene[0]["strain"] = "S"
    S_gene[0]["interaction"] = "monoculture"

    S_gene[1]["gene"] = gene
    S_gene[1]["strain"] = "S"
    S_gene[1]["interaction"] = "monoculture"
    
    # run monoculture E0 simulations
    E0_gene = alpha_table["E"][i]
    repressed_E0 = make_scaled_model(E0_iM, E0_gene, gene)
    E_gene = run_standard_monoculture_simulation(model = repressed_E0, carbon = "lcts_e", carbon_mmol = carbon_conc, auxotrophies = "met__L_e", auxo_mmol = met_conc, max_cycles = num_single_repression_cycles)

    # clean S monoculture dataframes
    E_gene[0]["gene"] = gene
    E_gene[0]["strain"] = "E"
    E_gene[0]["interaction"] = "monoculture"

    E_gene[1]["gene"] = gene
    E_gene[1]["strain"] = "E"
    E_gene[1]["interaction"] = "monoculture"
        
    # run coculture simulation
    ES_gene = run_standard_coculture_simulation(model1 = repressed_E0, model2 = repressed_S0, carbon = "lcts_e", carbon_mmol = carbon_conc, max_cycles = num_single_repression_cycles)

    ES_gene[0] = ES_gene[0].melt(
        id_vars=['cycle'],
        value_vars=[col for col in ES_gene[0].columns if col.startswith('iM') or col.startswith('STM')],
        var_name='strain',
        value_name='biomass'
    )

    ES_gene[0]["gene"] = gene
    ES_gene[0]["interaction"] = "coculture"

    ES_gene[1] = pivot_long_fluxes(ES_gene[1])
    ES_gene[1]["gene"] = gene
    ES_gene[1]["interaction"] = "coculture"
    ES_gene[1]["strain"] = "E"

    ES_gene[2] = pivot_long_fluxes(ES_gene[2])
    ES_gene[2]["gene"] = gene
    ES_gene[2]["interaction"] = "coculture"
    ES_gene[2]["strain"] = "S"

    # save individual data
    biomass = pd.concat([S_gene[0], E_gene[0], ES_gene[0]])
    fluxes = pd.concat([S_gene[1], E_gene[1], ES_gene[1], ES_gene[2]])
    exchanges = fluxes[fluxes['reaction'].str.contains('EX_', case=True, na=False) |
        fluxes['reaction'].str.contains('biomass', case=False, na=False)]

    # concatenate data
    biomass.to_csv(base + "biomasses/" + (gene + ".csv"))
    exchanges.to_csv(base + "fluxes/" + (gene + ".csv"))

    # print
    print("loop", i + 1, "of", len(gene_list) + 1)
    print("gene KO: ", gene)
