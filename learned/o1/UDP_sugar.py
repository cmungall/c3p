"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to a sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for uridine monophosphate (UMP) moiety
    ump_smarts = "n1(c(=O)[nH]c(=O)c1)[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O"  # Uridine with ribose
    ump_pattern = Chem.MolFromSmarts(ump_smarts)
    if ump_pattern is None:
        return False, "Unable to create uridine pattern"

    # Check for uridine moiety
    if not mol.HasSubstructMatch(ump_pattern):
        return False, "Uridine moiety not found"

    # Define SMARTS pattern for diphosphate linkage
    diphosphate_smarts = "OP(=O)(O)OP(=O)(O)O"  # Diphosphate group
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if diphosphate_pattern is None:
        return False, "Unable to create diphosphate pattern"

    # Check for diphosphate linkage connected to uridine
    ump_matches = mol.GetSubstructMatches(ump_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    found_uridine_diphosphate = False
    for ump_match in ump_matches:
        for dp_match in diphosphate_matches:
            # Get the atom indices of the ribose 5' oxygen and diphosphate oxygen
            ribose_O5_prime = ump_match[5]  # Adjust index based on pattern
            diphosphate_O = dp_match[0]  # First oxygen in diphosphate
            # Check if ribose O5' is connected to diphosphate
            bond = mol.GetBondBetweenAtoms(ribose_O5_prime, diphosphate_O)
            if bond:
                found_uridine_diphosphate = True
                break
        if found_uridine_diphosphate:
            break
    if not found_uridine_diphosphate:
        return False, "Diphosphate linkage not connected to uridine"

    # Define SMARTS pattern for sugar moiety (pyranose ring with multiple hydroxyls)
    sugar_smarts = "[C@H]1([O])[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"  # Simplified sugar ring
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return False, "Unable to create sugar pattern"

    # Check for sugar moiety
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Sugar moiety not found"

    # Check if sugar is connected to diphosphate
    found_sugar_linkage = False
    for dp_match in diphosphate_matches:
        for sugar_match in sugar_matches:
            # Get the atom indices of diphosphate oxygen and sugar anomeric carbon
            diphosphate_O = dp_match[-1]  # Last oxygen in diphosphate
            anomeric_C = sugar_match[0]  # First carbon in sugar ring
            # Check if diphosphate O is connected to sugar anomeric carbon
            bond = mol.GetBondBetweenAtoms(diphosphate_O, anomeric_C)
            if bond:
                found_sugar_linkage = True
                break
        if found_sugar_linkage:
            break
    if not found_sugar_linkage:
        return False, "Sugar moiety not connected via diphosphate linkage"

    return True, "UDP-sugar found with UDP moiety linked via diphosphate to sugar"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:46211',
        'name': 'UDP-sugar',
        'definition': 'A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.',
        'parents': []
    },
    'config': {}
}