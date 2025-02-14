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

    # Look for the oxime group pattern: C=N-O
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OX1H0]")
    if not mol.HasSubstructMatch(oxime_pattern):
        return False, "No oxime group found"

    # Check for aliphatic chain: ensure no aromatic carbons
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6):
        return False, "Contains aromatic carbons, not fully aliphatic"

    return True, "Contains an oxime group derived from an aliphatic aldehyde"

# Example test cases
example_smiles = [
    "OC(C(O)C(O)\\C=N\\O)C(O)CO", # aliphatic aldoxime
    "[H]\\C(=N/O)C(C)CC",         # aliphatic aldoxime
    "c1ccccc1C=NO",               # aromatic, not aliphatic
]

for smiles in example_smiles:
    result, reason = is_aliphatic_aldoxime(smiles)
    print(f"SMILES: {smiles}, Is Aliphatic Aldoxime: {result}, Reason: {reason}")