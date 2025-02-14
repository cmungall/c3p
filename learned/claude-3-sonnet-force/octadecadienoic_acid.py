"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:36334 octadecadienoic acid
Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly 18 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Found {c_count} carbon atoms, need exactly 18"

    # Check for exactly 2 double bonds
    n_double_bonds = rdMolDescriptors.CalcNumAromaticRings(mol) + rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_double_bonds != 2:
        return False, f"Found {n_double_bonds} double bonds, need exactly 2"

    # Check for straight chain
    chain_pattern = Chem.MolFromSmarts("[CH2,CH3]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Not a straight chain"

    return True, "Meets the criteria for an octadecadienoic acid"