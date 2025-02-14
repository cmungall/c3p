"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a molecule with a D-glucose moiety linked to another molecule.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
         bool: True if the molecule is a D-glucoside, False otherwise.
         str: Reason for the classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for D-glucose. Stereochemistry is specified.
    # This accounts for the ring form of D-glucose.
    d_glucose_smarts = '[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)'
    d_glucose_pattern = Chem.MolFromSmarts(d_glucose_smarts)

    #Check for a match to the D-glucose pattern.
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "No D-glucose moiety found"

    # Check for a glycosidic bond (an oxygen attached to the anomeric carbon of the glucose pattern)
    # Find the anomeric carbon (C1) in the glucose pattern
    glucose_match = mol.GetSubstructMatch(d_glucose_pattern)
    if glucose_match:
         anomeric_carbon = glucose_match[0]
    else:
        return False, "D-glucose pattern not found or invalid"
    
    glucose_atom = mol.GetAtomWithIdx(anomeric_carbon)

    has_glycosidic_bond = False
    for neighbor in glucose_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8: # Check for Oxygen
              #Check if it is connected to the glucose residue, not an internal oxygen
            is_glycosidic = False
            for neighbor_of_neighbor in neighbor.GetNeighbors():
                if neighbor_of_neighbor.GetIdx() not in glucose_match:
                    is_glycosidic = True
            if is_glycosidic:
               has_glycosidic_bond = True
               break

    if not has_glycosidic_bond:
      return False, "No glycosidic bond found"

    return True, "Contains a D-glucose moiety linked via a glycosidic bond."