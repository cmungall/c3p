"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is derived from an aliphatic aldehyde with an oxime group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect C=NOH oxime group pattern
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OH]")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No oxime group derived from an aldehyde (correct pattern applied)"

    # Exclude any aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Aromatic structures are present, thus not aliphatic"

    # Check for presence of aliphatic carbon chains
    # Assumes at least one decent-length aliphatic chain attached to C=NOH group
    aliphatic_pattern = Chem.MolFromSmarts("CC")
    if not mol.HasSubstructMatch(aliphatic_pattern):
        return False, "No sufficient aliphatic carbon chain found"

    # Ensure no ketone-derived oxime pattern, such as C=NOHR
    ketoxime_pattern = Chem.MolFromSmarts("[CX3](=[O])[R](=[NX2]-[OH])")
    if mol.HasSubstructMatch(ketoxime_pattern):
        return False, "Ketoxime pattern found, not an aldoxime"

    return True, "Contains an oxime group derived from an aliphatic aldehyde"

# Example test cases
example_smiles = [
    "OC(C(O)C(O)\\C=N\\O)C(O)CO", # (1E)-2,3,4,5,6-pentahydroxyhexanal oxime
    "[H]\\C(=N/O)C(C)CC",         # (E)-2-methylbutanal oxime
    "c1ccccc1C=NO",               # aromatic, not aliphatic
]

for smiles in example_smiles:
    result, reason = is_aliphatic_aldoxime(smiles)
    print(f"SMILES: {smiles}, Is Aliphatic Aldoxime: {result}, Reason: {reason}")