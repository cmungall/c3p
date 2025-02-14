"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition: Any organic compound having an initial boiling point less than or equal to 250°C measured at a standard atmospheric pressure of 101.3 kPa.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    Uses the Joback method to estimate the boiling point.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is organic (contains carbon and hydrogen)
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if 6 not in atom_nums:
        return False, "Molecule does not contain carbon"
    if 1 not in atom_nums:
        return False, "Molecule does not contain hydrogen"

    # Estimate boiling point using Joback method
    try:
        boiling_point = estimate_boiling_point_joback(mol)
    except Exception as e:
        return False, f"Boiling point estimation failed: {e}"

    if boiling_point <= 250:
        return True, f"Estimated boiling point is {boiling_point:.1f}°C which is ≤ 250°C"
    else:
        return False, f"Estimated boiling point is {boiling_point:.1f}°C which is > 250°C"

def estimate_boiling_point_joback(mol):
    """
    Estimates the boiling point of a molecule using the Joback method.

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object

    Returns:
        float: Estimated boiling point in °C
    """
    # Joback group contributions for boiling point
    joback_groups = {
        # Format: SMARTS pattern: (Delta Tb term, group description)
        "[CX4H3]":      (23.58, 'CH3'),
        "[CX4H2][!#1]": (22.88, 'CH2'),
        "[CX4H][!#1][!#1]": (21.74, 'CH'),
        "[CX4][!#1][!#1][!#1]": (20.23, 'C'),
        "[OX2H]":       (9.88, 'OH (alcohol)'),
        "[NX3H2]":      (16.54, 'NH2'),
        "[NX3H][!#1]":  (14.68, 'NH'),
        "[NX3][!#1][!#1]": (12.85, 'N'),
        "[$([CX3]=[OX1])][#6]": (7.73, 'C=O (ketone)'),
        "[$([CX3]=[OX1])][OX2H1]": (7.73, 'C=O (carboxylic acid)'),
        "[OX2][CX3](=O)[#6]": (7.73, 'C-O (ester)'),
        "[#6][OX2][#6]": (9.44, 'C-O-C (ether)'),
        "[#6][SX2][#6]": (9.44, 'C-S-C (thioether)'),
        "[SX2H]":       (9.44, 'S-H'),
        "[#6][F]":      (-10.44, 'C-F'),
        "[#6][Cl]":     (6.96, 'C-Cl'),
        "[#6][Br]":     (8.68, 'C-Br'),
        "[#6][I]":      (10.40, 'C-I'),
        # Add more group contributions as needed
    }

    # Sum of group contributions
    delta_tb = 198.0  # Base value for Tb calculation
    for smarts, (contribution, desc) in joback_groups.items():
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt)
        n_matches = len(matches)
        delta_tb += contribution * n_matches

    return delta_tb  # Boiling point in °C

# Example usage:
if __name__ == "__main__":
    smiles_list = [
        "CCCCCCCCC(O)CCCCC",  # octadecan-6-ol
        "ClC=C(Cl)Cl",        # trichloroethene
        "CCCCCCC",            # heptane
        "CC1CCCC1",           # methylcyclopentane
        "CCC(O)CC",           # pentan-3-ol
        "ClC(Cl)([H])[H]",    # dichloromethane-d2
    ]

    for smiles in smiles_list:
        result, reason = is_volatile_organic_compound(smiles)
        print(f"SMILES: {smiles}")
        print(f"Is VOC: {result}")
        print(f"Reason: {reason}")
        print("-" * 40)