"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a molecule composed of a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for recognition of structural components
    ceramide_pattern = Chem.MolFromSmarts("[C][C]([NH2])CC(=O)N[C@@H]")
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([O][C@@H]([C@H]([C@H]([C@@H]1O)O)O)CO)")
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@@H]1C[C@H]([NH2]C(=O)C)O[C@@H](C(O)=O)[C@H]1O")

    # Check for ceramide presence
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for sugar chains
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2:
        return False, f"Insufficient sugar moieties, found {len(sugar_matches)}"

    # Check for one or more sialic acids
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid found in the structure"

    # Check if the structure is fully connected as expected in gangliosides
    # Assuming simplified check to ensure linkage integrity
    if not all([mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in [ceramide_pattern, sugar_pattern, sialic_acid_pattern]]):
        return False, "Structure not recognized as a connected ganglioside"

    return True, "Contains ceramide backbone with attached glycosphingolipid structure and sialic acid"