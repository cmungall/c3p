"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:36357 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is defined as any molecule that contains four amino-acid residues 
    connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 4 peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 4:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 4"
    
    # Check for 4 amino acid residues
    # Use a generic pattern to match any alpha-amino acid (N-C(C)(C)-C(=O))
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][C@H]([C,c])([C,c])C(=O)"
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) != 4:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need exactly 4"
    
    # Check for N-terminus and C-terminus
    n_terminus_pattern = Chem.MolFromSmarts("[NX3H2]")
    c_terminus_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    n_terminus_match = mol.GetSubstructMatches(n_terminus_pattern)
    c_terminus_match = mol.GetSubstructMatches(c_terminus_pattern)
    if len(n_terminus_match) != 1 or len(c_terminus_match) != 1:
        return False, "Missing N-terminus or C-terminus"
    
    return True, "Molecule contains 4 amino acid residues connected by peptide bonds"