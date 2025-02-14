"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:51139 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid carrying a hydroxy group on the
    carbon atom at position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for dicarboxylic acid pattern (two -C(=O)O groups)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 2:
        return False, "Not a dicarboxylic acid"
    
    # Find alpha carbon (between two carboxylic acids)
    alpha_carbon_idx = None
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            neighbors1 = [n.GetAtomicNum() for n in atom1.GetNeighbors()]
            neighbors2 = [n.GetAtomicNum() for n in atom2.GetNeighbors()]
            if 8 in neighbors1 and 8 in neighbors2:
                alpha_carbon_idx = atom1.GetIdx() if 8 in neighbors1 else atom2.GetIdx()
                break
    
    if alpha_carbon_idx is None:
        return False, "No alpha carbon found between two carboxylic acids"
    
    # Check for hydroxy group on alpha carbon
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    has_hydroxyl = any(neighbor.GetAtomicNum() == 8 and
                       neighbor.GetTotalNumHs() == 1
                       for neighbor in alpha_carbon.GetNeighbors())
    
    if not has_hydroxyl:
        return False, "No hydroxyl group on alpha carbon"
    
    return True, "Dicarboxylic acid with hydroxyl group on alpha carbon"