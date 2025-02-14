"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is characterized by a long hydrocarbon chain containing at least
    one carbon-carbon double bond (C=C), with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group: [CX3](=O)[OX1H0-,OX2H1]
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check that the molecule is not a conjugate or complex derivative
    # Exclude molecules with phosphate, sulfate, or sugar moieties
    unwanted_groups = [
        Chem.MolFromSmarts("P(=O)(O)(O)"),  # Phosphate group
        Chem.MolFromSmarts("S(=O)(=O)O"),   # Sulfate group
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")  # Sugar ring
    ]
    for group in unwanted_groups:
        if mol.HasSubstructMatch(group):
            return False, "Molecule contains groups not typical of fatty acids (e.g., phosphate, sulfate, sugars)"

    # Check for the presence of at least one carbon-carbon double bond: [C]=[C]
    c_c_double_bond = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(c_c_double_bond):
        return False, "No carbon-carbon double bond found (no C=C)"

    # Count the total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Too few carbon atoms ({c_count}), requires at least 8"

    # Check that the molecule is predominantly hydrocarbon (allowing for O in functional groups)
    # Count heteroatoms other than oxygen
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6, 8)]
    if len(heteroatoms) > 1:
        return False, f"Contains {len(heteroatoms)} heteroatoms other than oxygen, which is atypical for fatty acids"

    # Ensure the majority of bonds are single or double bonds between carbons
    # Calculate the fraction of carbons in the molecule
    total_atoms = mol.GetNumAtoms()
    if c_count / total_atoms < 0.5:
        return False, "Molecule is not predominantly hydrocarbon"

    # Allow small rings (e.g., epoxides), but exclude larger rings and aromatic systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Check if the rings are small (3 members) and contain oxygen (e.g., epoxides)
        for ring in ring_info.AtomRings():
            if len(ring) > 3:
                return False, "Molecule contains rings larger than epoxides, not typical for fatty acids"
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetAtomicNum() == 6 for atom in ring_atoms):
                return False, "Molecule contains hydrocarbon rings, not typical for fatty acids"
            if any(atom.GetIsAromatic() for atom in ring_atoms):
                return False, "Molecule contains aromatic rings, not typical for fatty acids"

    # Check for ester or amide linkages that would indicate the molecule is a derivative
    ester = Chem.MolFromSmarts("C(=O)O[C,c]")
    amide = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(ester):
        return False, "Molecule is an ester derivative of a fatty acid"
    if mol.HasSubstructMatch(amide):
        return False, "Molecule is an amide derivative of a fatty acid"

    # Check for multiple carboxylic acid groups (dicarboxylic acids are less typical)
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(carboxylic_matches) > 1:
        return False, "Molecule has multiple carboxylic acid groups, not a typical fatty acid"

    return True, "Molecule is an olefinic fatty acid with at least one C=C double bond and appropriate structure"