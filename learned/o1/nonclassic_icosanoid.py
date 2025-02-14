"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signaling molecule made by oxygenation
    of C20 fatty acids other than the classic icosanoids (the leukotrienes and the prostanoids).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Molecule has {c_count} carbons, expected 20"

    # Get total number of oxygen atoms
    o_count_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check for terminal carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No terminal carboxylic acid group found"

    # Count oxygens in terminal carboxylic acid group(s)
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    o_count_carboxy = len(carboxy_matches) * 2  # Each carboxylic acid has 2 oxygens

    # Calculate number of oxygens beyond terminal carboxylic acid(s)
    o_count_additional = o_count_total - o_count_carboxy
    if o_count_additional < 1:
        return False, "No additional oxygen-containing functional groups found"

    # Exclude molecules with cyclopentane ring (prostaglandins)
    cyclopentane = Chem.MolFromSmarts('C1CCCC1')
    if mol.HasSubstructMatch(cyclopentane):
        return False, "Molecule contains a cyclopentane ring (possible prostaglandin)"

    # Exclude molecules with conjugated triene system (leukotrienes)
    conjugated_triene = Chem.MolFromSmarts('C=C-C=C-C=C')
    triene_matches = mol.GetSubstructMatches(conjugated_triene)
    if triene_matches:
        # Further analysis could be added here to distinguish leukotrienes
        return False, "Molecule contains conjugated triene system (possible leukotriene)"

    # Check for oxygenated functional groups
    # Hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Epoxides (three-membered cyclic ethers)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    num_epoxides = len(mol.GetSubstructMatches(epoxide_pattern))

    # Ketones (C=O not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
    num_ketones = len(mol.GetSubstructMatches(ketone_pattern))

    # Check if there is at least one oxygenated functional group
    num_oxygenated_groups = num_hydroxyls + num_epoxides + num_ketones
    if num_oxygenated_groups < 1:
        return False, "No oxygenated functional groups beyond carboxylic acid found"

    # If molecule passes all checks, classify as nonclassic icosanoid
    return True, "Molecule meets criteria for nonclassic icosanoid"