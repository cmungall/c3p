"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide contains four amino-acid residues connected by peptide linkages.

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

    # Define a more flexible SMARTS pattern for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[!#1]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3 for tetrapeptide"

    # SMARTS pattern to detect amino acid residues (any N-C-C(=O) fragment with diversity in side chains)
    aa_residue_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)")
    aa_residue_matches = mol.GetSubstructMatches(aa_residue_pattern)
    
    if len(aa_residue_matches) != 4:
        return False, f"Detected {len(aa_residue_matches)} amino acid residues, need exactly 4"
    
    # Verify connection: ensuring four amino acids are connected sequentially
    # This assumes the identified peptide bonds and amino acids residues define a linear, fully connected molecule.
    connected_components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(connected_components) != 1:
        return False, "Molecule not fully connected as a single chain"
    
    return True, "Contains four amino-acid residues connected by three peptide bonds, valid tetrapeptide structure"