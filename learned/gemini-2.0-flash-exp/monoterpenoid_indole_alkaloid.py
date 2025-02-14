"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety, a monoterpenoid unit and is an alkaloid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for Indole Substructure:
    indole_pattern = Chem.MolFromSmarts("c1c[nH]c2ccccc12")
    if not mol.HasSubstructMatch(indole_pattern):
         return False, "No indole substructure found"

    # 2. Check for Complexity:
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Too few rings for a monoterpenoid indole alkaloid."

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbons for a monoterpenoid indole alkaloid."

    # Count nitrogens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "Must have at least one nitrogen (alkaloid nature)"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight is too low for a monoterpenoid indole alkaloid."
    

    return True, "Contains an indole substructure, multiple rings, and sufficient complexity"