"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:36689 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is defined as a tetraterpenoid (C40) with a polyene chain (alternating C=C and C-C bonds)
    and may contain cyclic end groups and/or oxygen atoms (xanthophylls).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyene chain (minimum 8 conjugated double bonds)
    polyene_pattern = Chem.MolFromSmarts("[C]=[C][C]=[C][C]=[C][C]=[C][C]=[C][C]=[C][C]=[C]")
    polyene_matches = mol.GetSubstructMatches(polyene_pattern)
    if not polyene_matches:
        return False, "No polyene chain found (minimum 8 conjugated double bonds)"

    # Check for cyclic end groups
    cyclic_end_pattern = Chem.MolFromSmarts("[C@]1(CCCC1)[C]=[C]")
    cyclic_end_matches = mol.GetSubstructMatches(cyclic_end_pattern)

    # Check for oxygen atoms (xanthophylls)
    oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check for typical carotenoid backbone (C40 or close)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, "Carbon count outside typical carotenoid range (35-45)"

    # Classify based on findings
    if cyclic_end_matches and oxygens:
        return True, "Contains polyene chain, cyclic end groups, and oxygens (xanthophyll carotenoid)"
    elif cyclic_end_matches:
        return True, "Contains polyene chain and cyclic end groups (carotene carotenoid)"
    elif oxygens:
        return True, "Contains polyene chain and oxygens (xanthophyll carotenoid)"
    else:
        return True, "Contains polyene chain (carotene carotenoid)"