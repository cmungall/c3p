"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise.
        str: Reason for classification.
    """

    # Handle the simple proton case first
    if smiles == "[H+]":
         return True, "Single proton cation"

    #split into molecules, handle each separately
    smiles_list = smiles.split(".")
    
    is_overall_cation = False
    reason = "Not a cation"

    for s in smiles_list:

        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return False, "Invalid SMILES string"

        net_charge = 0
        for atom in mol.GetAtoms():
            net_charge += atom.GetFormalCharge()

        if net_charge > 0:
              is_overall_cation = True
              reason = "Contains at least one positive charge and no negative charges."
        elif net_charge < 0:
             return False, "Not a cation, net negative charge."
    
    return is_overall_cation, reason