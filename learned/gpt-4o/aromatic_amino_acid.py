"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid includes an amino acid structure with an aromatic ring attached as a side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern to identify amino acids: N and C connected with further chain continuity
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=O)[OX2H]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid functional group found"

    # Aromatic ring pattern
    aromatic_ring_pattern = Chem.MolFromSmarts("a")
    
    # We look for aromatic rings attached to the alpha carbon or adjacent to it
    aromatic_side_chain_pattern = Chem.MolFromSmarts("[NX3][CX4H]([CX3](=O)[OX2H])[CX4,CX3](a)")
    
    # Check if there is any aromatic side chain pattern attached to amino acids 
    if mol.HasSubstructMatch(aromatic_side_chain_pattern):
        return True, "Contains amino acid functional group with an aromatic ring side chain"
    
    return False, "No aromatic ring attached as a side chain to an amino acid"