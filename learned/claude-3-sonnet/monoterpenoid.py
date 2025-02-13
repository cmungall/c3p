"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:24064 monoterpenoid

A monoterpenoid is any terpenoid derived from a monoterpene. The term includes
compounds in which the C10 skeleton of the parent monoterpene has been rearranged
or modified by the removal of one or more skeletal atoms (generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of atoms
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10 or n_atoms > 20:
        return False, "Atom count outside the range for monoterpenoids"

    # Check for C10 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Less than 10 carbon atoms"

    # Check for rearranged or modified monoterpene skeleton
    monoterpene_skeleton_pattern = Chem.MolFromSmarts("[C&r6,C&r5]")
    if mol.HasSubstructMatch(monoterpene_skeleton_pattern):
        return True, "Contains rearranged or modified monoterpene skeleton"

    # Check for characteristic functional groups and substructures
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    ether_pattern = Chem.MolFromSmarts("COC")
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(ketone_pattern) or mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains characteristic functional group or substructure of monoterpenoids"

    # Check for long hydrocarbon chains (but not as a strict rule)
    hydrocarbon_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    if mol.HasSubstructMatch(hydrocarbon_pattern):
        return False, "Likely a sesquiterpenoid or diterpenoid"

    # Default case: not enough evidence to classify
    return None, None