import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
import itertools



# This SMARTS does the carboligation of the two alpha-keto acids - it's a really long string
def rule_1():
    rule1 = '[#6:1](=[#8:2])[#6:3](=[#8:4])[#8:5].[#6:6](=[#8:7])[#6:8](=[#8:9])[#8:10]>>[O:2]=[C:1][C:6]([O:7])[C:8](=[O:9])[O:10]'
    return rule1



# This SMARTS does the conversion from the intermediate to the final product
def rule_2():
    rule2 = '[O:2]=[C:1][C:6]([O:7])[C:8](=[O:9])[O:10]>>[C:1][C:6]'
    return rule2



# Making this rule 3, for in the event that carbon dioxide falls off of intermediate 1
# Keith said that this could happen
def rule_3():
    #fomerly had it output CO2 as well [O:1]=[C;H0:2]([O;H1:3])[C;H0:4][O;H1:5]>>[O:1]=[C:2]=[O:3].[C:4][O:5]
    rule3 = '[O:1]=[C;H0:2]([O;H1:3])[C;H0:4][O;H1:5]>>[C:4][O:5]'
    return rule3


def make_mol(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return molecule



def make_smiles(mol):
    smiles = Chem.MolToSmiles(mol)
    return smiles



"""rule must be in SMARTS form, molecules must be SMILES"""
def run_reaction_1(reactant1, reactant2, rule):
    # takes in reactants 1&2 in SMILES form and converts them to mol objects
    mol_reactant_1 = make_mol(reactant1)
    mol_reactant_2 = make_mol(reactant2)

    # then it turns both of these mol objects into a tuple, as this is what 'RunReactants requires'
    react_tuple_1 = (mol_reactant_1, mol_reactant_2)
    react_tuple_2 = (mol_reactant_2, mol_reactant_1)
    reaction = Chem.ReactionFromSmarts(rule)

    output_1 = reaction.RunReactants(react_tuple_1)
    output_2 = reaction.RunReactants(react_tuple_2)
    output_list = [output_1,output_2]

    return_list = []

    for each_output in output_list:
        for each_rxn in each_output:
            for each_product in each_rxn:
                return_list.append(Chem.MolToSmiles(each_product))

    return return_list



def run_reaction_2(reactant, rule):
    # takes in reactants 1&2 in SMILES form and converts them to mol objects
    mol_reactant_1 = make_mol(reactant)

    # then it turns both of these mol objects into a tuple, as this is what 'RunReactants requires'
    react_tuple = (mol_reactant_1,)
    reaction = Chem.ReactionFromSmarts(rule)

    output = reaction.RunReactants(react_tuple)

    for each_rxn in output:
        for each_product in each_rxn:
            return Chem.MolToSmiles(each_product)


"""Will give you back the intermediate, but without a CO2, which may fall off during rxn 1"""
def run_reaction_3(reactant,rule):
    # takes in reactants 1&2 in SMILES form and converts them to mol objects
    mol_reactant_1 = make_mol(reactant)

    # then it turns both of these mol objects into a tuple, as this is what 'RunReactants requires'
    react_tuple = (mol_reactant_1,)
    reaction = Chem.ReactionFromSmarts(rule)

    output = reaction.RunReactants(react_tuple)

    for each_rxn in output:
        for each_product in each_rxn:
            return Chem.MolToSmiles(each_product)



def get_monoisotopic_mass(smiles):
    return Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles))



def output_file(dataframe):
    """this is my function
    :param dataframe this is my dataframe
    :type dataframe pds
    :return nothing
    :rtype NAN
    """
    filepath = '/Users/grad/Documents/IBiS/Tyo_Lab/Alpha_Keto_Acids/Alpha_Keto_Output.csv'
    dataframe.to_csv(filepath)



def main():
    df = pd.read_csv('/Users/grad/Documents/IBiS/Tyo_Lab/Alpha_Keto_Acids/Alpha_Keto_Inputs.csv')
    df2 = pd.DataFrame(columns=['Reactant_1_SMILES',
                                'Reactant_1_Name',
                                'Reactant_1_Monoisotopic_Mass',
                                'Reactant_2_SMILES',
                                'Reactant_2_Names',
                                'Reactant_2_Monoisotopic_Mass',
                                'Intermediates',
                                'Intermediates_Monoisotopic_Mass',
                                'Final_Product',
                                'CO2_Fall_Off_Intermediate',
                                'CO2_Fall_Off_Intermediate_Monoisotopic_Mass'])
    df2.index.name = 'Index'

    rule1 = rule_1()
    rule2 = rule_2()
    rule3 = rule_3()

    iterlist = []
    for i in range(df.shape[0]):
        iterlist.append(i)

    counter = 0
    iterate = itertools.combinations_with_replacement(iterlist, 2)
    for num in iterate:
        name_1 = df.iloc[num[0],0]
        name_2 = df.iloc[num[1],0]
        alpha_keto_1 = df.iloc[num[0], 1]
        alpha_keto_2 = df.iloc[num[1], 1]

        #getting monoisotopic masses for alpha_keto_1 & alpha_keto_2
        alpha_keto_1_monoisotopic_mass = get_monoisotopic_mass(alpha_keto_1)
        alpha_keto_2_monoisotopic_mass = get_monoisotopic_mass(alpha_keto_2)

        reaction_1_output = run_reaction_1(alpha_keto_1, alpha_keto_2, rule1)
        intermediate1 = reaction_1_output[0]
        intermediate2 = reaction_1_output[1]
        #intermediates = [intermediate1,intermediate2]
        temp_intermediates = [intermediate1,intermediate2]
        temp_intermediates2 = set(temp_intermediates)
        intermediates = ';'.join(temp_intermediates2)

        monoisotopic_mass1 = get_monoisotopic_mass(intermediate1)
        #monoisotopic_mass2 = get_monoisotopic_mass(intermediate2)

        final_product = run_reaction_2(intermediate1, rule2)

        # adding this in in case a CO2 falls of the intermediate
        fall_off_intermediate_1 = run_reaction_3(intermediate1,rule3)
        fall_off_intermediate_2 = run_reaction_3(intermediate2,rule3)
        temp_fall_off = [fall_off_intermediate_1,fall_off_intermediate_2]
        temp_fall_off_2 = set(temp_fall_off)
        fall_offs = ';'.join(temp_fall_off_2)

        #get the fall offs monoisotopic mass for intermediate 1
        fall_offs_monoisotopic_mass = get_monoisotopic_mass(fall_off_intermediate_1)


        append_me = [alpha_keto_1,
                     name_1,
                     round(alpha_keto_1_monoisotopic_mass,4),
                     alpha_keto_2,
                     name_2,
                     round(alpha_keto_2_monoisotopic_mass,4),
                     intermediates,
                     round(monoisotopic_mass1,4),
                     final_product,
                     fall_offs,
                     round(fall_offs_monoisotopic_mass,4)]

        df2.loc[counter] = append_me
        counter += 1

    output_file(df2)



if __name__ == '__main__':
    main()