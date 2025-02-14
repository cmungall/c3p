"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36363 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid scaffold
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@H]3[C@@H]([C@@H]1CC[C@@H]2[C@@H]4[C@H]([C@H]3[C@@H](O)CC4)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid scaffold found"

    # Check for 3-hydroxy group in beta position
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H]1CCC2(C)CC[C@H]([C@@H]3[C@H]([C@H]2[C@@H](C1)C)CCC4=CC(=O)CC[C@]34C)C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "3-hydroxy group not in beta position"

    # Check for typical steroid properties
    num_rings = rdMolDescriptors.CalcNumRingAtoms(mol)
    if num_rings < 4:
        return False, "Not enough rings for a steroid"

    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 1:
        return False, "Too many aromatic rings for a steroid"

    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 200 or mol_weight > 600:
        return False, "Molecular weight outside typical range for steroids"

    return True, "Contains steroid scaffold with 3-hydroxy group in beta position"