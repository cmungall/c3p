"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.
    It should not be a peptide (no peptide bonds).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    carboxylic_acid_smarts = "[CX3](=O)[OX1H0-,OX2H1]"  # Carboxylic acid (include deprotonated form)
    amino_group_smarts = "[NX3;!$(N=*);!$([N+])]"  # Any nitrogen not double-bonded or positively charged
    peptide_bond_smarts = "[CX3](=O)[NX3][CX4]"  # Peptide bond: O=C-N-C (amide bond connected to carbon)
    amide_bond_smarts = "[CX3](=O)[NX3]"  # General amide bond

    # Create pattern molecules
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    amino_group_pattern = Chem.MolFromSmarts(amino_group_smarts)
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    amide_bond_pattern = Chem.MolFromSmarts(amide_bond_smarts)
    
    # Search for carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    else:
        num_carboxylic_acids = len(carboxylic_acid_matches)
    
    # Search for amino groups
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if not amino_group_matches:
        return False, "No amino group found"
    else:
        num_amino_groups = len(amino_group_matches)
    
    # Identify nitrogens involved in amide bonds
    amide_nitrogens = set()
    amide_matches = mol.GetSubstructMatches(amide_bond_pattern)
    for match in amide_matches:
        amide_nitrogens.add(match[1])  # Index of nitrogen in amide bond
    
    # Identify free amino groups (not involved in amide bonds)
    free_amino_groups = []
    for match in amino_group_matches:
        nitrogen_idx = match[0]
        if nitrogen_idx not in amide_nitrogens:
            free_amino_groups.append(nitrogen_idx)
    
    if not free_amino_groups:
        return False, "No free amino group found (all involved in amide bonds)"
    else:
        num_free_amino_groups = len(free_amino_groups)
    
    # Identify carbons involved in peptide bonds
    peptide_carbons = set()
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    for match in peptide_matches:
        peptide_carbons.add(match[0])  # Index of carbonyl carbon in peptide bond
    
    # Identify free carboxylic acid groups (not involved in peptide bonds)
    free_carboxylic_acids = []
    for match in carboxylic_acid_matches:
        carbon_idx = match[0]
        if carbon_idx not in peptide_carbons:
            free_carboxylic_acids.append(carbon_idx)
    
    if not free_carboxylic_acids:
        return False, "No free carboxylic acid group found (all involved in peptide bonds)"
    else:
        num_free_carboxylic_acids = len(free_carboxylic_acids)
    
    # Optionally, exclude large molecules (e.g., with many amide bonds)
    total_amide_bonds = len(amide_matches)
    if total_amide_bonds > 1:
        return False, "Contains multiple amide bonds, likely a peptide"
    
    return True, f"Contains {num_free_carboxylic_acids} free carboxylic acid group(s) and {num_free_amino_groups} free amino group(s)"