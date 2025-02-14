"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has a terminal aldehyde group and a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal aldehyde group (C=O, connected to one other C)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) != 1:
        return False, f"Found {len(aldehyde_matches)} terminal aldehyde groups, need exactly 1"

    # Check for long carbon chain (at least 6 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
       return False, "Carbon chain too short"

    # Check for number of carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for fatty aldehyde"
    if o_count != 1:
        return False, "Must have exactly 1 oxygen from the aldehyde group"

    
    # Check molecular weight is greater than 80. A shorter chain may have MW close to 80,
    # therefore more weight is given to presence of at least 6 carbons.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80:
         return False, "Molecular weight too low for fatty aldehyde"
        
    
    return True, "Has a terminal aldehyde group and a long carbon chain"