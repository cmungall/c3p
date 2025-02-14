"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:16670 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule containing two amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for exactly two peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 2"
    
    # Check for at least one amino acid side chain
    side_chain_pattern = Chem.MolFromSmarts("[C&!$(C=O)][NX3H2]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if not side_chain_matches:
        return False, "No amino acid side chains found"
    
    # Check for carboxyl and amino groups
    has_carboxyl = any(atom.GetSymbol() == "O" and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in mol.GetAtomNeighbors(atom.GetIdx())) > 0 for atom in mol.GetAtoms())
    has_amino = any(atom.GetSymbol() == "N" and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in mol.GetAtomNeighbors(atom.GetIdx())) > 1 for atom in mol.GetAtoms())
    if not (has_carboxyl and has_amino):
        return False, "Missing carboxyl or amino group"
    
    return True, "Contains two amino acid residues connected by peptide linkages"