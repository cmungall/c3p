"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to find a peptide bond
    #  N-C(=O)-C where N has at least one H, to avoid other amides. 
    peptide_bond_pattern = Chem.MolFromSmarts("[N;!H0]-C(=O)-[#6]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    num_peptide_bonds = len(peptide_bonds)
    
    # Check for the minimum number of peptide bonds
    # A polypeptide with 10 amino acid residues will have at least 9 peptide bonds.
    if num_peptide_bonds < 9:
        return False, f"Molecule contains {num_peptide_bonds} peptide bonds, which is less than the required 9 for a polypeptide with 10 or more amino acid residues."
    
    return True, f"Molecule contains {num_peptide_bonds} peptide bonds, consistent with a polypeptide."