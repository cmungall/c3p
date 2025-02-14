"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:17855 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphatidylinositol backbone
    pi_pattern = Chem.MolFromSmarts("C(OC[C@@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)OC)=O")
    if not mol.HasSubstructMatch(pi_pattern):
        return False, "No phosphatidylinositol backbone found"
    
    # Check for at least one phosphate group on the inositol ring
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1OP(O)(O)=O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No phosphate groups on the inositol ring"
    
    # Count the number of phosphate groups attached to the inositol ring
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    n_phosphates = 0
    for match in phosphate_matches:
        atom = mol.GetAtomWithIdx(match[0])
        if atom.GetNeighbors()[0].GetNeighbors()[0].IsInRingSize(6):
            n_phosphates += 1
    
    if n_phosphates < 1:
        return False, "No phosphate groups attached to the inositol ring"
    
    return True, f"Contains phosphatidylinositol backbone with {n_phosphates} phosphate group(s) on the inositol ring"