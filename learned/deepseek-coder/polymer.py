"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a mixture composed of macromolecules of different kinds, which may be differentiated by composition, length, degree of branching, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for large molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a polymer"

    # Count the number of rotatable bonds to assess chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds for a polymer"

    # Check for a minimum number of atoms
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 50:
        return False, "Too few atoms for a polymer"

    # Look for specific repeating patterns common in polymers
    repeating_patterns = [
        Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]"),  # Carbon chain
        Chem.MolFromSmarts("[OX2]~[CX4]~[OX2]"),  # Ether linkage
        Chem.MolFromSmarts("[CX3](=[OX1])~[OX2]~[CX4]"),  # Ester linkage
        Chem.MolFromSmarts("[NX3]~[CX4]~[NX3]"),  # Amide linkage
        Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]"),  # Longer carbon chain
        Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]"),  # Even longer carbon chain
    ]
    
    repeating_matches = 0
    for pattern in repeating_patterns:
        if mol.HasSubstructMatch(pattern):
            repeating_matches += 1

    if repeating_matches < 2:
        return False, "No repeating patterns found"

    # If all checks pass, classify as polymer
    return True, "Contains repeating units, high molecular weight, and complex structure indicative of a polymer"