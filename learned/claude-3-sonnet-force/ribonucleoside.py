"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:27638 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for nucleobase
    nuc_base_pattern = Chem.MolFromSmarts("*nc*")
    if not mol.HasSubstructMatch(nuc_base_pattern):
        return False, "No nucleobase found"
    
    # Check for ribose sugar
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H](CO)O[C@H](O)[C@@H]1O)O")
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No ribose sugar found"
    
    # Check for N-glycosidic bond between base and ribose
    atom_idx_pairs = mol.GetSubstructMatches(Chem.MolFromSmarts("[C@H]1[C@@H](CO)O[C@H](O)[C@@H]1O.n"))
    for ribose_idx, base_idx in atom_idx_pairs:
        ribose_atom = mol.GetAtomWithIdx(ribose_idx)
        base_atom = mol.GetAtomWithIdx(base_idx)
        if mol.GetBondBetweenAtoms(ribose_idx, base_idx):
            return True, "Contains a nucleobase and ribose sugar connected via N-glycosidic bond"
    
    return False, "No N-glycosidic bond found between nucleobase and ribose"