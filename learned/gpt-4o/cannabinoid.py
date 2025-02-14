"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are identified by characteristic core structures and functional group profiles.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expand and refine patterns for known cannabinoid core structures
    cannabinoid_core_patterns = [
        Chem.MolFromSmarts('c1cc(Cc2ccccc2)ccc1'),  # Biphenyl structures
        Chem.MolFromSmarts('[C@H]1CCC(C)=C[C@H]1'),  # Cyclohexene with stereochemistry
        Chem.MolFromSmarts('Cc1cc(O)c(O)cc1'),       # Phenol with additional OH or methyl groups
    ]
    
    has_core_structure = any(mol.HasSubstructMatch(pattern) for pattern in cannabinoid_core_patterns)
    if has_core_structure:
        # Refine check with characteristic functional groups part of cannabinoids
        characteristic_groups = [
            Chem.MolFromSmarts('[OX2H]'),             # Hydroxyl group
            Chem.MolFromSmarts('O-[CX3](=O)'),        # Ester (carboxylate)
            Chem.MolFromSmarts('NC(=O)'),             # Amide
            Chem.MolFromSmarts('C=C(CC)CC'),          # Long alkenyl or alkyl chains
        ]
        
        has_characteristic_group = any(mol.HasSubstructMatch(fg) for fg in characteristic_groups)
        if has_characteristic_group:
            return True, "Characteristic cannabinoid patterns detected"
    
    return False, "No characteristic patterns of cannabinoids detected"