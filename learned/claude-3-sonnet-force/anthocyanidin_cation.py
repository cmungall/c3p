"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:36334 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an oxygenated derivative of flavylium (2-phenylchromenylium) with a positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavylium core (chromen-4-ylium)
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c2ccccc2cc1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core found"
    
    # Look for phenyl ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"
    
    # Look for oxygens outside the flavylium core
    oxygen_pattern = Chem.MolFromSmarts("O")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    oxygen_counts = [0] * mol.GetNumAtoms()
    for match in oxygen_matches:
        oxygen_counts[match] += 1
    n_oxygens_outside = sum(1 for count in oxygen_counts if count > 0 and not mol.GetAtomWithIdx(i).IsInRingSize(6))
    if n_oxygens_outside < 1:
        return False, "No oxygens outside flavylium core"
    
    # Check for positive charge
    formal_charge = rdMolDescriptors.CalcFormalCharge(mol)
    if formal_charge != 1:
        return False, "Incorrect formal charge (expected +1)"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15 or c_count > 25:
        return False, "Carbon count outside typical range (15-25)"
    if o_count < 4 or o_count > 10:
        return False, "Oxygen count outside typical range (4-10)"
    
    return True, "Contains flavylium core with phenyl ring and oxygens outside the core, with +1 charge"