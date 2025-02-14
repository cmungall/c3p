"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component
    attached to an unspecified sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): A tuple containing:
            - True if molecule is a UDP-sugar, False otherwise.
            - Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify the UDP core
    # Define a SMARTS pattern for the UDP core. This includes the uracil, ribose and diphosphate.
    udp_core_pattern = Chem.MolFromSmarts('O=C1NC(=O)C=CN1[C@H]2[C@@H](O)[C@@H](O)[C@H](COP(=O)(O)OP(=O)(O)O)[C@H]2O')
    if not mol.HasSubstructMatch(udp_core_pattern):
       return False, "No UDP core found"


    # 2. Verify the diphosphate linkage
    # Check that the terminal phosphate of UDP is connected to a sugar moiety.
    # We'll check for a P-O-C bond where the C is not part of the UDP core
    linkage_pattern = Chem.MolFromSmarts("[#15](=[OX1])(-[OX2])-[CX4]")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    if not linkage_matches:
       return False, "No diphosphate linkage found"

    # Check if the carbon involved in the linkage is part of the sugar by searching for a pattern
    # of [CH2,CH][O] at least 3 times
    sugar_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4][OX2][CX4][OX2]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
      return False, "No sugar moiety found"

    # 3. Check the overall composition
    # Check the correct number of atoms of each type
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if n_count < 2:
      return False, "Too few nitrogens"
    if o_count < 10:
        return False, "Too few oxygens"
    if c_count < 10:
      return False, "Too few carbons"
    if p_count != 2:
      return False, "Must have exactly 2 phosphorus atoms"

    # If all checks passed, it's likely a UDP-sugar
    return True, "Contains UDP core with diphosphate linkage to a sugar"