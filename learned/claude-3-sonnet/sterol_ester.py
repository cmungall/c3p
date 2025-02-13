"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C][C]2[C]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check if the ester is at position 3 of the steroid core
    # The pattern looks for a steroid core with an ester at position 3
    sterol_ester_pattern = Chem.MolFromSmarts("[C]1[C][C]([OX2][CX3](=[OX1]))[C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C][C]2[C]1")
    if not mol.HasSubstructMatch(sterol_ester_pattern):
        return False, "Ester group not at position 3 of steroid core"

    # Additional checks for typical sterol characteristics
    
    # Count carbons (sterols typically have >20 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol ester"
    
    # Look for typical sterol side chain (at least 2 carbons)
    side_chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH](C)*")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "Missing typical sterol side chain"
    
    # Count rings (sterols typically have 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for a sterol structure"
        
    return True, "Contains steroid core with ester group at position 3"