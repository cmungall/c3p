"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a glycerol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for two phosphatidic acid patterns, characterized by a phosphate group connected to carbon chains
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O[C@H]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"
    
    # Confirm the presence of long carbon chains (indicative of fatty acids) for each phosphatidic acid
    carbon_chain_pattern = Chem.MolFromSmarts("C(C(=O)O)C")  # A simple pattern to identify esters in the chains
    ester_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(ester_matches) < 4:
        return False, f"Insufficient ester groups for two phosphatidic acids, found {len(ester_matches)}"

    # Verify a central glycerol backbone connecting the phosphates (C-C-C with oxygens)
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O)C(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No central glycerol structure found"

    return True, "Molecule structure is consistent with cardiolipin"