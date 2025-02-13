"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:24702 polymer
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdqueries
from typing import Tuple

def is_polymer(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a mixture, which is composed of macromolecules of different kinds
    and which may be differentiated by composition, length, degree of branching, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of repeating structural units
    pattern = Chem.MolFromSmarts("*~*~*~*")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No repeating structural units found"

    # Check for high molecular weight (>500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for polymer"

    # Check for high number of rotatable bonds (>10)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient flexibility for polymer"

    # Check for presence of common polymer functional groups
    polymer_smarts = ['[NX3+]', '[OX2H]', '[OX1H0-]', '[CX3](=[OX1])[OX2H]', '[SX4+2]', '[PX4]']
    polymer_groups_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in polymer_smarts)
    if not polymer_groups_found:
        return False, "No common polymer functional groups found"

    # Check for high degree of branching
    branching_count = sum(1 for atom in mol.GetAtoms() if len(atom.GetNeighbors()) >= 3)
    if branching_count < 5:
        return False, "Low degree of branching for polymer"

    return True, "Molecule exhibits properties of a polymer"