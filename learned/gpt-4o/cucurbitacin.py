"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Looking for cucurbitane backbone structure
    tetracyclic_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCCC4')
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "No tetracyclic cucurbitane backbone structure found"

    # Checking for characteristic functional groups, e.g., hydroxyl (OH) and keto groups (C=O)
    hydroxyl_count = len(re.findall(r'\[OH\]', smiles))
    keto_pattern = Chem.MolFromSmarts('[CX3](=O)')
    keto_matches = mol.GetSubstructMatches(keto_pattern)

    # Cucurbitacins often possess multiple hydroxyl groups and keto moieties (at least 2 keto groups)
    if hydroxyl_count < 2 or len(keto_matches) < 2:
        return False, f"Insufficient characteristic functional groups: hydroxyls {hydroxyl_count}, ketones {len(keto_matches)}"

    # Check the molecular weight - cucurbitacins typically have higher molecular weights
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a typical cucurbitacin"

    return True, "Contains cucurbitane-like tetracyclic pattern with characteristic functional groups"