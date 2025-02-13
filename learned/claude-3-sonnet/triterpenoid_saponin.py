"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:36367 triterpenoid saponin
A terpene glycoside in which the terpene moiety is a triterpenoid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for triterpenoid backbone (6 isoprene units)
    triterpenoid_pattern = Chem.MolFromSmarts("[C@@]1([C@H](C)[C@H]2[C@@]3(C)CC[C@](O)([C@@H]3CC[C@@]12C)C)C")
    if not mol.HasSubstructMatch(triterpenoid_pattern):
        return False, "No triterpenoid backbone found"
    
    # Look for glycoside (sugar) groups attached
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]([C@H](O)[C@H](O)[C@H](O)O)O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar groups found"
    
    # Check molecular weight - triterpenoid saponins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for triterpenoid saponin"
    
    # Count carbons, oxygens, and hydrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if c_count < 30:
        return False, "Too few carbons for triterpenoid saponin"
    if o_count < 5:
        return False, "Too few oxygens for triterpenoid saponin"
    if h_count < 40:
        return False, "Too few hydrogens for triterpenoid saponin"
    
    return True, "Contains a triterpenoid backbone with attached sugar groups"