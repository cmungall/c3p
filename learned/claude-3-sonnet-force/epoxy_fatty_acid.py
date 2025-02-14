"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:46787 epoxy fatty acid
A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide pattern (O1C2CCCCC2C1)
    epoxide_pattern = Chem.MolFromSmarts("[O;R]1[C@;@R]2[C@@;@R][C@@;@R][C@;@R][C@@;@R][C@;@R]2[C@;@R]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Look for carboxylic acid pattern (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for long carbon chains (fatty acid)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be fatty acids"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 10:
        return False, "Too few carbons for fatty acid"
    if o_count < 2:
        return False, "Too few oxygens for epoxy fatty acid"

    return True, "Contains an epoxide ring and a carboxylic acid group with long fatty acid chains"