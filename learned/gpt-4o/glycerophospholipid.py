"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is defined as having a glycerol backbone with a phosphate group 
    ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphate group pattern (P=O(OH)2 or variants with ester linkages)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Ensure the phosphate is linked to a terminal carbon of the glycerol backbone
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    phosphate_atoms = mol.GetSubstructMatch(phosphate_pattern)

    for glycerol_atom in [glycerol_atoms[0], glycerol_atoms[2]]:  # Check both terminal carbons
        for phosphate_atom in phosphate_atoms:
            if mol.GetBondBetweenAtoms(glycerol_atom, phosphate_atom) is not None:
                # We found a glycerol-phosphate link
                return True, "Contains glycerol backbone with phosphate group attached"
    
    return False, "Phosphate group not ester-linked to glycerol terminal carbon"