"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16beta-hydroxy steroid
A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for stereochemistry
    mol = Chem.AddHs(mol)
    
    # Basic steroid core pattern (four fused rings)
    # Using SMARTS that matches the basic steroid skeleton
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"
    
    # Pattern for 16-beta-hydroxy group
    # The [H] specifies explicit hydrogen, OH is the hydroxy group
    # The @ symbols specify the stereochemistry
    beta_oh_pattern = Chem.MolFromSmarts("[C]12[C][C@H](O)[CH2][C@]1([CH2,CH3])[C@@H]3[C][C][C]2")
    
    if not mol.HasSubstructMatch(beta_oh_pattern):
        return False, "No 16-beta-hydroxy group found"
    
    # Additional checks for reasonable molecular weight and atom counts
    # Steroids typically have molecular weights between 250-1000
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 1000:
        return False, f"Molecular weight {mol_weight} outside typical steroid range (250-1000)"
    
    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:
        return False, f"Too few carbons ({carbon_count}) for a steroid structure"
    
    # Check for at least one oxygen (for the hydroxy group)
    oxygen_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if oxygen_count < 1:
        return False, "No oxygen atoms found"
        
    return True, "Contains steroid core with 16-beta-hydroxy group"