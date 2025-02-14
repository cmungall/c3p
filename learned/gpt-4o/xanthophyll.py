"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid characterized by a long chain
    of conjugated double bonds and one or more oxygen-containing functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for extensive conjugated double bonds system
    conjugated_pattern = Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]")  # C=C-C=C-C=C
    conjugated_systems = mol.GetSubstructMatches(conjugated_pattern)
    
    if len(conjugated_systems) < 3:  # A reasonable threshold, e.g., 3 matches
        return False, "Insufficient conjugated system typical for carotenoids"
    
    # Check for oxygen-containing functional groups
    oxy_functional_groups = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl
    oxy_groups = mol.GetSubstructMatches(oxy_functional_groups)
    
    if len(oxy_groups) < 1:
        return False, "No significant oxygen-containing groups found"

    # Extra checks for functionality diversity:
    # Keto group
    keto_functional_groups = Chem.MolFromSmarts("[CX3](=O)[#6]")
    keto_groups = mol.GetSubstructMatches(keto_functional_groups)

    # Epoxy group
    epoxy_functional_groups = Chem.MolFromSmarts("[OX2R]([#6])[#6]")
    epoxy_groups = mol.GetSubstructMatches(epoxy_functional_groups)
    
    if len(oxy_groups) + len(keto_groups) + len(epoxy_groups) < 2:
        return False, "Lacks diverse oxygen functional groups"

    # Check for at least 9 double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 9:
        return False, f"Insufficient number of conjugated double bonds; found {double_bond_count}"
    
    return True, "Chemical structure fits the typical xanthophyll characteristics of conjugated double bonds and diverse oxygen functionality"