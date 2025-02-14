"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two free carboxy groups (-COOH or -COO^-),
    not involved in ester or amide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxyl group pattern (both protonated and deprotonated forms)
    carboxyl_group_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    if carboxyl_group_pattern is None:
        return False, "Failed to create carboxyl group pattern"

    # Find all carboxyl groups in the molecule
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_group_pattern)
    num_carboxyl_groups = 0

    for match in carboxyl_matches:
        carboxyl_carbon_idx = match[0]
        hydroxyl_oxygen_idx = match[1]

        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        hydroxyl_oxygen = mol.GetAtomWithIdx(hydroxyl_oxygen_idx)

        # Check if hydroxyl oxygen is connected only to carboxyl carbon and possibly hydrogen (not ester)
        oxygen_neighbors = [a.GetAtomicNum() for a in hydroxyl_oxygen.GetNeighbors() if a.GetIdx() != carboxyl_carbon_idx]
        if any(atomic_num > 1 for atomic_num in oxygen_neighbors):
            continue  # It's part of an ester

        # Check if carboxyl carbon is not single-bonded to nitrogen (not amide)
        carbon_neighbors = [a for a in carboxyl_carbon.GetNeighbors() if a.GetIdx() != hydroxyl_oxygen_idx]
        for neighbor in carbon_neighbors:
            if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                bond = mol.GetBondBetweenAtoms(carboxyl_carbon_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    break  # It's part of an amide
        else:
            # Increment count if not part of ester or amide
            num_carboxyl_groups += 1

    # Exclude peptides by checking for peptide bonds (amide bonds between carbonyl carbon and nitrogen)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Molecule contains peptide bonds"

    # Exclude large molecules (e.g., molecular weight > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight is too high ({mol_wt:.2f} Da)"

    # Check if there are exactly two carboxyl groups
    if num_carboxyl_groups == 2:
        return True, "Molecule contains exactly two free carboxyl groups"
    else:
        return False, f"Found {num_carboxyl_groups} free carboxyl group(s), expected exactly 2"