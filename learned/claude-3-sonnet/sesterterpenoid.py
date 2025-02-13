"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:36236 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are terpenoids derived from a sesterterpene (C25 skeleton) which may
    have been rearranged or modified by removing skeletal atoms (e.g., methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of atoms
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 25:
        return False, "Less than 25 atoms, too small for sesterterpenoid"
    
    # Check for characteristic C25 skeleton
    c25_pattern = Chem.MolFromSmarts("[C&r5]")  # Ring of size 5 containing only carbons
    c25_matches = mol.GetSubstructMatches(c25_pattern)
    if not c25_matches:
        return False, "No C25 skeleton found"
    
    # Look for terpenoid features
    has_isoprene_unit = mol.HasSubstructMatch(Chem.MolFromSmarts("C=CC(C)C"))
    if not has_isoprene_unit:
        return False, "No isoprene units found, not a terpenoid"
    
    # Check for modifications to parent skeleton
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return True, "Contains C25 skeleton with few rotatable bonds, likely sesterterpenoid"
    else:
        return True, "Contains C25 skeleton and isoprene units, with modifications - sesterterpenoid"