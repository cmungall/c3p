"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:38756 N-acetyl-amino acid
An N-acyl-amino acid that has acetyl as the acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetyl group attached to nitrogen
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group attached to nitrogen"
    
    # Look for amino acid backbone (C-N-C-C=O)
    amino_acid_pattern = Chem.MolFromSmarts("[C;X4][N;X3][C;X4][C;X3](=[O;X1])[O;X2]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"
    
    # Check if the nitrogen attached to the acetyl group is also part of the amino acid backbone
    acetyl_N = mol.GetSubstructMatches(acetyl_pattern)[0][1]
    amino_acid_N = mol.GetSubstructMatches(amino_acid_pattern)[0][1]
    if acetyl_N != amino_acid_N:
        return False, "Acetyl group not attached to amino acid nitrogen"
    
    # Count number of nitrogens - should be 1 for N-acetyl-amino acids
    n_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogens != 1:
        return False, "More than one nitrogen atom present"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    return True, "Contains acetyl group attached to nitrogen of an amino acid backbone"