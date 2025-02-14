"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of four amino acid residues connected by three peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of the amino acid residue backbone pattern
    amino_acid_pattern = Chem.MolFromSmarts("N[CX4][CX3](=O)")
    residue_matches = mol.GetSubstructMatches(amino_acid_pattern)

    if len(residue_matches) < 4: #must have at least four
         return False, f"Found {len(residue_matches)} amino acid residues, need at least 4"

    #check for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=[O])-[N]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 3:
      return False, f"Found {len(peptide_bonds)} peptide bonds, need at least 3"

    return True, "Contains four amino acid residues connected by at least three peptide linkages"