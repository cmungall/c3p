"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:36851 nonclassic icosanoid

A nonclassic icosanoid is any biologically active signalling molecule made by oxygenation
of C20 fatty acids other than the classic icosanoids (the leukotrienes and the prostanoids).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, "Does not have 20 carbon atoms"

    # Check for oxygenation (at least 2 oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Not sufficiently oxygenated"

    # Check for a carbon chain of at least 10 carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No long carbon chain found"

    # Check for at least 2 double bonds (unsaturated)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Not sufficiently unsaturated"

    # Check for at least 2 oxygenated functional groups (hydroxyl, epoxide, hydroperoxide, ketone)
    oxygenated_pattern = Chem.MolFromSmarts("[OX1H,OX2,OX3,OX4]")
    oxygenated_matches = mol.GetSubstructMatches(oxygenated_pattern)
    if len(oxygenated_matches) < 2:
        return False, "Not enough oxygenated functional groups"

    # Check for absence of prostanoid and leukotriene patterns
    prostanoid_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@@H](C1)C(=O)C[C@H]2")
    leukotriene_pattern = Chem.MolFromSmarts("[C@@H]1[C@H](C[C@H]1C(=O)O)C(=O)O")
    if mol.HasSubstructMatch(prostanoid_pattern) or mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Molecule matches classic icosanoid pattern"

    return True, "Molecule is a C20 fatty acid derivative with oxygenation and unsaturation"