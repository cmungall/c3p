"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a glucuronic acid moiety with a glycosidic bond at the C1 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucuronic acid substructure with specific stereochemistry
    # using SMARTS. We are checking the correct orientation of the hydroxyl groups
    # and the carboxyl group, and the ring oxygen
    glucuronic_acid_smarts = "[C@H]1([O])[C@@H]([O])[C@@H]([O])[C@H]([C@@H]1O)C(=O)O"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid substructure found"

    # Get the matches for the substructure
    matches = mol.GetSubstructMatches(glucuronic_acid_pattern)

    #Check that one of the C1 is attached to an oxygen
    glycosidic_bond_found = False
    for match in matches:
       c1_index = match[0] # the index of C1 in the molecule
       c1_atom = mol.GetAtomWithIdx(c1_index)
       for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8: # check for the glycosidic bond with Oxygen
                glycosidic_bond_found = True
                break
       if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
       return False, "No glycosidic bond at C1 of glucuronic acid."


    return True, "Contains beta-D-glucuronic acid with a glycosidic bond at the C1 position."