"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35654 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    Secondary alpha-hydroxy ketones have a carbonyl group and a hydroxy group linked
    to the same carbon, which is also connected to an organyl group and a hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the alpha-hydroxy ketone pattern: C(=O)C(O)
    pattern = Chem.MolFromSmarts("[C&D2&H1]")
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No secondary alpha-hydroxy ketone group found"
    
    # Make sure the carbon is not part of a ring
    for atom_idx in matches:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.IsInRing():
            return False, "Alpha-hydroxy ketone group is part of a ring"
    
    # Check for exactly 1 carbonyl and 1 hydroxyl group
    carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2)
    
    if carbonyl_count != 1 or hydroxyl_count != 1:
        return False, "Incorrect number of carbonyl/hydroxyl groups"
    
    # Check molecular weight - secondary alpha-hydroxy ketones typically <500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for secondary alpha-hydroxy ketone"
    
    return True, "Contains a secondary alpha-hydroxy ketone group"