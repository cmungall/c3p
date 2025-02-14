"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:18020 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins have two methyl substituents at position 2, a carboxylate substituent
    at position 3, and a carboxamido group at position 6 of the penam ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for penam ring system (4-membered beta-lactam ring fused to a 5-membered ring)
    penam_pattern = Chem.MolFromSmarts("[C&R1]=4[C&R1]=3[C&R2]=2[C&R2]=1[N&R3]=4[C&R1]=3[C&R4]=2[C&R1]=1")
    if not mol.HasSubstructMatch(penam_pattern):
        return False, "No penam ring system found"
    
    # Check for two methyl groups at position 2
    methyl_pattern = Chem.MolFromSmarts("[C&R2]([C])(C)")
    if len(mol.GetSubstructMatches(methyl_pattern)) != 1:
        return False, "Incorrect number of methyl groups at position 2"
    
    # Check for carboxylate at position 3
    carboxylate_pattern = Chem.MolFromSmarts("[C&R3](=O)[O-]")
    if len(mol.GetSubstructMatches(carboxylate_pattern)) != 1:
        return False, "No carboxylate group at position 3"
    
    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("[N&R4][C](=O)[O]")
    if len(mol.GetSubstructMatches(carboxamido_pattern)) != 1:
        return False, "No carboxamido group at position 6"
    
    # Additional checks (not strictly required)
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 15 or n_atoms > 50:
        return False, "Molecule size atypical for penicillins"
    
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings != 2:
        return False, "Incorrect number of rings for penicillins"
    
    return True, "Contains penam ring system with two methyl groups at position 2, carboxylate at position 3, and carboxamido group at position 6"