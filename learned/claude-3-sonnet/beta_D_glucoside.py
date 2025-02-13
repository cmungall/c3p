"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:36688 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is any D-glucoside where the anomeric center has a beta configuration.

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
    
    # Look for D-glucose substructure
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose substructure found"
    
    # Check for beta configuration at anomeric center
    anomeric_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if any(n.GetSymbol() == "C" and n.GetDegreeTotal() == 4 for n in neighbors):
                anomeric_atom = atom
                break
    
    if anomeric_atom is None:
        return False, "Could not find anomeric center"
    
    # Check configuration at anomeric center
    config = mol.GetAtomWithIdx(anomeric_atom.GetIdx()).GetProp("_CIPCode")
    if config != "S":
        return False, "Anomeric center is not beta configuration"
    
    return True, "Molecule contains a beta-configured D-glucose moiety"