"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:27475 Icosanoid
Icosanoids are signalling molecules arising from oxidation of the three C20 essential fatty acids (EFAs) 
icosapentaenoic acid (EPA), arachidonic acid (AA) and dihomo-gamma-linolenic acid (DGLA).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Found {c_count} carbon atoms, icosanoids must have exactly 20"

    # Look for at least one cyclopentane ring (characteristic of prostaglandins)
    cyclopentane_pattern = Chem.MolFromSmarts("[C@H]1CCC[C@H]1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        cyclopentene_pattern = Chem.MolFromSmarts("[C@H]1=CCC[C@H]1")
        if not mol.HasSubstructMatch(cyclopentene_pattern):
            return False, "No cyclopentane or cyclopentene ring found (characteristic of prostaglandins)"

    # Look for at least 3 carbon-carbon double bonds (characteristic of polyunsaturated fatty acids)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found only {len(double_bond_matches)} carbon-carbon double bonds, icosanoids should have at least 3"

    # Look for at least one oxygen (characteristic of oxidized fatty acids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, icosanoids should be oxidized fatty acids"

    # Check for long carbon chains (characteristic of fatty acids)
    longest_chain = max([len(chain) for chain in AllChem.GetMolFrags(mol)])
    if longest_chain < 10:
        return False, "Longest carbon chain too short to be a fatty acid"

    return True, "Contains characteristics of an icosanoid: 20 carbons, cyclopentane/ene ring, multiple C=C bonds, oxygens, and long carbon chains"