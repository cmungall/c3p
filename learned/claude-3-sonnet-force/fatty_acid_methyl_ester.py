"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:35473 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a carboxylic ester obtained by the formal condensation
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for methyl ester group (-C(=O)-O-C)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"

    # Look for long carbon chain (at least 6 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Carbon chain too short for fatty acid"

    # Check molecular formula for expected elements (C, H, O)
    allowed_atoms = set([6, 1, 8])  # C, H, O
    atom_nums = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if not atom_nums.issubset(allowed_atoms):
        return False, "Unexpected atoms present, should only contain C, H, O"

    return True, "Contains a methyl ester group and a long carbon chain (fatty acid)"