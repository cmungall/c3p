"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33709 alpha-amino acid

An alpha-amino acid is defined as an amino acid in which the amino group is 
located on the carbon atom at the position alpha to the carboxy group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino group and carboxyl group on adjacent carbons
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1]")
    
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    for amino_idx in amino_matches:
        amino_atom = mol.GetAtomWithIdx(amino_idx)
        amino_neighbors = [n.GetIdx() for n in amino_atom.GetNeighbors()]
        
        for carboxyl_idx in carboxyl_matches:
            carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
            carboxyl_neighbors = [n.GetIdx() for n in carboxyl_atom.GetNeighbors()]
            
            # Check if amino and carboxyl groups are on adjacent carbons
            common_neighbor = set(amino_neighbors) & set(carboxyl_neighbors)
            if common_neighbor:
                alpha_carbon_idx = common_neighbor.pop()
                alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
                if alpha_carbon.GetAtomicNum() == 6:  # Carbon
                    return True, "Contains amino group and carboxyl group on adjacent carbons (alpha-amino acid)"
    
    return False, "No alpha-amino acid substructure found"