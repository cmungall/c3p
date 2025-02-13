"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: CHEBI:51274 long-chain fatty alcohol

A fatty alcohol with a chain length ranging from C13 to C22.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"
    
    # Check for long carbon chain (C13-C22)
    carbon_chain_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "Carbon chain shorter than C13"
    
    carbon_chain_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if carbon_chain_matches:
        return False, "Carbon chain longer than C22"
    
    # Count atoms and check for saturation
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if (h_count - 1) != (2*c_count + 2):  # -1 for OH, +2 for terminal CH3
        return False, "Unsaturated hydrocarbon chain"
    
    return True, "Molecule is a long-chain fatty alcohol (C13-C22)"