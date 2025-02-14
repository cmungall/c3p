"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:24192 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is any D-glucoside in which the anomeric center has beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for D-glucose substructure
    d_glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H](C([C@@H]([C@@H]1O)O)O)O)O)O")
    d_glucose_matches = mol.GetSubstructMatches(d_glucose_pattern)
    if not d_glucose_matches:
        return False, "No D-glucose substructure found"

    # Check for beta-configuration at anomeric center
    for match in d_glucose_matches:
        anomeric_center = match[-1]  # Last atom in the match is the anomeric center
        atom = mol.GetAtomWithIdx(anomeric_center)
        if atom.GetHybridization() == Chem.HybridizationType.SP3:
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 8:  # Oxygen neighbor
                    cfg = mol.GetAtomProp(str(neighbor.GetIdx()), "_CIPCode")
                    if cfg == "R":  # beta-configuration
                        return True, "Contains D-glucose with beta-configuration at anomeric center"
    
    return False, "D-glucose substructure found but not in beta-configuration"