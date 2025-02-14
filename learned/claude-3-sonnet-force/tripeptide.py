"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:36357 tripeptide
Any oligopeptide that consists of three amino-acid residues connected by peptide linkages.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for 3 peptide bonds (-C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[C&x3](=O)[N&x3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"
    
    # Count amino acid residues (fragments connected by peptide bonds)
    amino_acid_pattern = Chem.MolFromSmarts("[N&x3]-[C&x3](=O)-[C&x3](=[O&x1])-[C&x3]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) != 3:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need exactly 3"
    
    # Check for N-terminus and C-terminus
    if mol.GetNumAtoms() < 6:
        return False, "Molecule too small to be a peptide"
    
    n_terminus = any(atom.GetSymbol() == "N" and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 2 for atom in mol.GetAtoms())
    c_terminus = any(atom.GetSymbol() == "C" and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    
    if not n_terminus or not c_terminus:
        return False, "Missing N-terminus or C-terminus"
    
    # Check for common protecting groups
    protecting_groups = ["Boc", "Cbz", "Fmoc", "Ac", "Bz", "Trt", "Z"]
    has_protecting_group = any(group in Chem.MolToSmiles(mol) for group in protecting_groups)
    
    if has_protecting_group:
        reason = "Contains 3 amino acid residues connected by peptide bonds with protecting group"
    else:
        reason = "Contains 3 amino acid residues connected by peptide bonds"
    
    return True, reason