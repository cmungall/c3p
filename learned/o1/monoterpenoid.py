"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is any terpenoid derived from a monoterpene, including compounds
    in which the C10 skeleton of the parent monoterpene has been rearranged or modified
    by the removal of one or more skeletal atoms (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define monoterpene core structures (common monoterpenes)
    monoterpene_cores = [
        Chem.MolFromSmarts("C=C(C)C[C@@H](C)CC=C(C)C"),  # Myrcene
        Chem.MolFromSmarts("C1=CC=C(C=C1)C(C)C"),         # p-Cymene
        Chem.MolFromSmarts("CC(=CCCC(=CC)C)C"),           # Ocimene
        Chem.MolFromSmarts("C1CC2CCC1C2C"),               # Bornane
        Chem.MolFromSmarts("CC1=CCC(CC1)C(C)C"),          # Limonene
        Chem.MolFromSmarts("CC1=CCC2C1CCC2C"),            # Pinene
        Chem.MolFromSmarts("CC(C)C1CCC(=CC1)C"),          # Sabinene
        Chem.MolFromSmarts("CC1CCC(CC1)C(C)(C)O"),        # Terpinen-4-ol
    ]

    # Check if the molecule contains any of the monoterpene core structures
    has_core = False
    for core in monoterpene_cores:
        if mol.HasSubstructMatch(core):
            has_core = True
            break

    if not has_core:
        return False, "No monoterpene core structure detected"

    # Monoterpenoids may contain various functional groups
    # Check for common modifications (functional groups)
    functional_group_patterns = {
        "alcohol": "[OX2H]",                # Hydroxyl group
        "ketone": "[CX3](=O)[#6]",          # Ketone group
        "aldehyde": "[CX3H1](=O)[#6]",      # Aldehyde group
        "carboxylic_acid": "C(=O)[OH1]",    # Carboxylic acid group
        "ester": "C(=O)O[#6]",              # Ester group
        "ether": "[OD2]([#6])[#6]",         # Ether group
        "epoxide": "[OX2r3]",               # Epoxide ring
        "peroxide": "OO",                   # Peroxide group
    }

    has_functional_group = False
    for fg_name, smarts in functional_group_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_functional_group = True
            break  # We can stop after finding one functional group

    if not has_functional_group:
        return False, "No common monoterpenoid functional groups detected"

    # Monoterpenoids may be cyclic or acyclic and may have rearranged skeletons
    # Additional checks can include ring structures and backbone modifications
    # For this implementation, we'll accept the molecule if it passes the above checks

    # Compile reasons for classification
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    num_rings = mol.GetRingInfo().NumRings()

    reason = f"Contains monoterpene core structure, {num_carbons} carbons, "
    reason += f"{num_oxygens} oxygens, "
    reason += f"{'cyclic' if num_rings > 0 else 'acyclic'} structure, "
    reason += "with common monoterpenoid functional groups"

    return True, reason