"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:37434 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is a ceramide consisting of sphinganine (a long-chain aminoalcohol)
    in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphinganine backbone pattern: NH-C-C-OH
    sphinganine_pattern = Chem.MolFromSmarts("[NH1][CX4][CX4][OH1]")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Look for amide bond (N-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"
    
    # Check for fatty acid chain (long carbon chain attached to amide)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chain too short to be fatty acid"
    
    # Count carbon atoms in fatty acid chain
    atom_ids = list(set([atom.GetIdx() for match in fatty_acid_matches for atom in match]))
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() in atom_ids)
    if c_count < 12:
        return False, "Fatty acid chain too short"

    return True, "Contains sphinganine backbone with fatty acyl group attached via amide bond"