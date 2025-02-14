"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a short chain of amino acids linked by peptide bonds.

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
    
    # Look for peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_matches) < 1: #An oligopeptide should have at least one peptide bond
        return False, f"Found {len(peptide_matches)} peptide bonds, at least 1 required"
    
    # Look for amino acid residue pattern (N-C-C(=O))
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    if len(amino_acid_matches) < 2 or len(amino_acid_matches) > 25: #set limits for oligopeptides
        return False, f"Found {len(amino_acid_matches)} amino acid residues, must be between 2 and 25"

    # Check for free carboxylic acid and amine groups (or amide bond)
    has_terminal_carboxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O"))
    has_terminal_amino    = mol.HasSubstructMatch(Chem.MolFromSmarts("N"))
    has_amide             = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3][CX3](=[OX1])"))

    if not (has_terminal_carboxyl or has_amide) or not (has_terminal_amino or has_amide):
        return False, "Not a typical oligopeptide structure without terminal amino or carboxyl group and without amides"
    
    # Check chain length using molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 3000: # set an upper limit on MW based on typical oligopeptide, increased
        return False, "Molecular weight is too high for a typical oligopeptide"

    #check for the presence of nitrogens and oxygens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if n_count < 1 or o_count < 2:
        return False, "Oligopeptide must have at least one N and at least two O atoms"
    
    return True, "Has peptide bonds and amino acid residues, within size range of a typical oligopeptide"