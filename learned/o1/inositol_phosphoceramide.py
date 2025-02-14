"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide consists of an inositol residue linked via a phosphodiester bridge
    to a ceramide moiety (sphingoid base linked via an amide bond to a fatty acyl chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define inositol ring pattern (cyclohexane ring with hydroxyl groups on each carbon)
    inositol_pattern = Chem.MolFromSmarts("C1(CO)C(CO)C(CO)C(CO)C(CO)O1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Define ceramide moiety pattern
    # Amide bond linking fatty acyl chain to sphingoid base (long-chain amino alcohol)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@H](CO)[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # Define phosphodiester bridge pattern
    phosphodiester_pattern = Chem.MolFromSmarts("P(=O)(O[C])[O][C]")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester bridge found"

    # Check connectivity between inositol ring and phosphate group
    # Get indices of matched atoms for inositol and phosphate
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    phosphate_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    inositol_atoms = set(sum(inositol_matches, ()))
    phosphate_atoms = set(sum(phosphate_matches, ()))

    # Check if inositol ring is connected to phosphate group
    connected = False
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in inositol_atoms and end_idx in phosphate_atoms) or \
           (begin_idx in phosphate_atoms and end_idx in inositol_atoms):
            connected = True
            break
    if not connected:
        return False, "Phosphate group is not connected to inositol ring"

    # Check connectivity between ceramide moiety and phosphate group
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    ceramide_atoms = set(sum(ceramide_matches, ()))

    connected = False
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in ceramide_atoms and end_idx in phosphate_atoms) or \
           (begin_idx in phosphate_atoms and end_idx in ceramide_atoms):
            connected = True
            break
    if not connected:
        return False, "Phosphate group is not connected to ceramide moiety"

    return True, "Contains inositol phosphoceramide structure"