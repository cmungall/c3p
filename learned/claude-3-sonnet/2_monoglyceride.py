"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with a single fatty acid chain 
    attached via an ester bond at position 2, and free hydroxyls at positions 1 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and standardize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens and kekulize
    mol = Chem.AddHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Pattern for terminal primary alcohols (more permissive)
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]")
    oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    if len(oh_matches) != 2:
        return False, f"Found {len(oh_matches)} terminal hydroxyl groups, need exactly 2"

    # Check for 2-position attachment using a simpler pattern
    # Matches central carbon with ester attachment
    position2_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4]([OX2][CX3]=[OX1])-[CH2X4]")
    if not mol.HasSubstructMatch(position2_pattern):
        return False, "Ester group not at position 2"

    # Verify presence of fatty acid chain (at least 4 carbons)
    # More flexible pattern that handles both saturated and unsaturated chains
    fatty_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No fatty acid chain found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 7:  # minimum 3 for glycerol + at least 4 for fatty acid
        return False, "Too few carbons for 2-monoglyceride"
    if o_count != 4:  # exactly 4 oxygens: 2 hydroxyls + 2 for ester group
        return False, "Must have exactly 4 oxygens (2 hydroxyls + 1 ester group)"

    return True, "Contains glycerol backbone with single fatty acid chain at position 2"