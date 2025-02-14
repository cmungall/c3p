"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments, with a benzopyran-4-one core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core benzopyran-4-one structure (flavone or flavonol)
    # Flavone (C=O at position 4):
    flavone_core_pattern = Chem.MolFromSmarts("c1ccc(-c2cc(=O)c3ccccc3o2)cc1") # more precise
    # Flavonol (C=O at position 4, and OH at position 3):
    flavonol_core_pattern = Chem.MolFromSmarts("c1ccc(-c2c(O)cc(=O)c3ccccc3o2)cc1")

    if not (mol.HasSubstructMatch(flavone_core_pattern) or mol.HasSubstructMatch(flavonol_core_pattern)):
        return False, "No benzopyran-4-one core structure found"
    
    # 2. Check for Hydroxyl and Methoxy groups. Count how many. At least one of each
    num_hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    num_methoxys = len(mol.GetSubstructMatches(Chem.MolFromSmarts("OC")))

    if num_hydroxyls == 0 and num_methoxys == 0:
        return False, "No hydroxyl or methoxy groups found"

    # 3. Check for glycosylation (O-C-O outside of the core pattern)
    glycoside_pattern = Chem.MolFromSmarts("O-C-O")
    core_match = mol.GetSubstructMatch(flavone_core_pattern) if mol.HasSubstructMatch(flavone_core_pattern) else mol.GetSubstructMatch(flavonol_core_pattern)
    if core_match:
      core_atoms = set(core_match)
      glycoside_matches = [match for match in mol.GetSubstructMatches(glycoside_pattern) if not set(match) <= core_atoms]
      if len(glycoside_matches) > 0:
         return True, "Benzopyran-4-one core structure with hydroxyl and/or methoxy, and possible glycosylation found"

    # 4. Check if additional aromatic rings are present
    aromatic_rings = Chem.MolFromSmarts("c1ccccc1")
    num_aromatic_rings = len(mol.GetSubstructMatches(aromatic_rings))

    if num_aromatic_rings >= 2:
        return True, "Benzopyran-4-one core structure with hydroxyl and/or methoxy and additional aromatic ring found"
    
    return True, "Benzopyran-4-one core structure with hydroxyl and/or methoxy found"