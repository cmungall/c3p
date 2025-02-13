"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:36357 tripeptide
Any oligopeptide that consists of three amino-acid residues connected by peptide linkages.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

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
    
    # Look for 2 amide bonds (-C(=O)-N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 2:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 2"
    
    # Look for 3 amino acid residues (N-C-C)
    aa_pattern = Chem.MolFromSmarts("NC(C)C")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    if len(aa_matches) != 3:
        return False, f"Found {len(aa_matches)} amino acid residues, need exactly 3"
    
    # Check for peptide bonds connecting residues
    peptide_bond_pattern = Chem.MolFromSmarts("N(C(=O))C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 2"
    
    # Check hydrogen bond donors/acceptors
    n_hbond_donors = rdMolDescriptors.CalcNumHBD(mol)
    n_hbond_acceptors = rdMolDescriptors.CalcNumHBA(mol)
    if n_hbond_donors < 3 or n_hbond_donors > 6:
        return False, "Hydrogen bond donor count outside expected range"
    if n_hbond_acceptors < 4 or n_hbond_acceptors > 8:
        return False, "Hydrogen bond acceptor count outside expected range"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside expected range for tripeptide"
    
    # Additional checks for structural diversity
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[r]")):  # Ring system
        return True, "Contains 3 amino acid residues linked by peptide bonds (cyclic)"
    if any(atom.GetDegree() > 4 for atom in mol.GetAtoms()):  # Quaternary atoms
        return True, "Contains 3 amino acid residues linked by peptide bonds (quaternary atoms present)"
    
    return True, "Contains 3 amino acid residues linked by peptide bonds"