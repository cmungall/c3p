"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:36322 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count amino acid residues
    pattern = Chem.MolFromSmarts("[N;X3;H2,H1&!$(N(-C=O)-O-C=O)][C;X4][C;X3](=O)[N;X3]")
    matches = mol.GetSubstructMatches(pattern)
    num_residues = len(matches)
    
    # Oligopeptides typically have 3-20 amino acid residues
    if num_residues < 3:
        return False, f"Only {num_residues} amino acid residues found, need at least 3"
    elif num_residues > 20:
        return False, f"{num_residues} amino acid residues found, too many for an oligopeptide"
    
    # Look for peptide bonds and amino acid sidechains
    peptide_bond_pattern = Chem.MolFromSmarts("[N;X3][C;X3](=[O;X1])[C;X4]")
    has_peptide_bonds = mol.HasSubstructMatch(peptide_bond_pattern)
    
    sidechain_pattern = Chem.MolFromSmarts("[N;X3][C;X4][C;X3](=O)[N;X3]")
    has_sidechains = mol.HasSubstructMatch(sidechain_pattern)
    
    if not has_peptide_bonds or not has_sidechains:
        return False, "Missing peptide bonds or amino acid sidechains"
    
    # Check molecular weight - oligopeptides typically 300-2000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for oligopeptide"
    elif mol_wt > 2000:
        return False, "Molecular weight too high for oligopeptide"
    
    return True, "Contains a small number of amino acid residues linked by peptide bonds"