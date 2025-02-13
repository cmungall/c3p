"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: CHEBI:36196 unsaturated fatty acid
Any fatty acid containing at least one C=C or C#C bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is a fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for alkene/alkyne bonds (C=C or C#C)
    unsaturated_pattern = Chem.MolFromSmarts("C=,#C")
    if not mol.HasSubstructMatch(unsaturated_pattern):
        return False, "No unsaturated bonds (C=C or C#C) found"
    
    # Count rotatable bonds to check for long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Molecule is too rigid, not a long fatty acid chain"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for fatty acid"
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (carboxylic acid)"
    
    return True, "Contains carboxylic acid group and at least one C=C or C#C unsaturation"