"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: CHEBI:36343 Cyclic fatty acid
A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify non-aromatic aliphatic rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    aliphatic_rings = [ring for ring in rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() == 0 for idx in ring)]

    if not aliphatic_rings:
        return False, "No non-aromatic aliphatic rings found in the structure"

    # Check for carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons for a fatty acid"

    # Additional checks or filters can be added here if needed

    return True, "Contains a carboxylic acid group and at least one non-aromatic aliphatic ring in the structure"