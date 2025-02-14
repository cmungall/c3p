"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is defined as an acyclic branched or unbranched hydrocarbon having
    exactly one carbon-carbon double bond with the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure only contains Carbon and Hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1):
            return False, "Molecule contains elements other than carbon and hydrogen"

    # Ensure exactly one carbon-carbon double bond (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    num_double_bonds = len(double_bond_matches)

    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} carbon-carbon double bonds, need exactly 1"

    # Use RDKit to calculate the molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    
    # Extract counts of C and H from the formula
    c_count = 0
    h_count = 0
    import re
    for match in re.finditer(r'([A-Z][a-z]*)(\d*)', formula):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        if element == 'C':
            c_count = count
        elif element == 'H':
            h_count = count
    
    if c_count == 0 or h_count == 0 or h_count != 2 * c_count:
        return False, f"Molecular formula does not match CnH2n (found C{c_count}H{h_count})"
    
    return True, "Molecule is an alkene with one carbon-carbon double bond and fits the general formula CnH2n"