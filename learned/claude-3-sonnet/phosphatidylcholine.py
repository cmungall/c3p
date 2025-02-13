"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: phosphatidylcholine
A glycerophosphocholine that is glycero-3-phosphocholine bearing two acyl substituents at positions 1 and 2.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphocholine group (-P(=O)([O-])OCC[N+](C)(C)C)
    phosphocholine_pattern = Chem.MolFromSmarts("[P](=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for phosphate group with negative charge
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No negatively charged phosphate group found"

    # Check for quaternary ammonium (choline) with positive charge
    choline_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
        return False, "No positively charged choline group found"

    # Check that the molecule has the correct number of charged groups
    pos_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    if pos_charges != 1 or neg_charges != 1:
        return False, f"Incorrect number of charges (need 1 positive, 1 negative; found {pos_charges} positive, {neg_charges} negative)"

    # Verify carbon chains (fatty acids) are present
    carbon_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = len(mol.GetSubstructMatches(carbon_chain))
    if chain_matches < 2:
        return False, "Missing fatty acid chains"

    return True, "Contains glycerol backbone with two fatty acid chains and phosphocholine group"