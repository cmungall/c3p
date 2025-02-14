"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl is characterized by a 1,4-diphenylbenzene skeleton and may have substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core para-terphenyl structure using SMARTS
    para_terphenyl_smarts = "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"
    para_terphenyl_pattern = Chem.MolFromSmarts(para_terphenyl_smarts)

    # Check if the structure contains the para-terphenyl skeleton
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains para-terphenyl skeleton (1,4-diphenylbenzene)"
    else:
        return False, "No para-terphenyl skeleton found"

# Example test
smiles_example = "O(C1=C(O)C(C2=CC=CC=C2)=CC(=C1C3=CC=CC=C3)O)C"
result, reason = is_para_terphenyl(smiles_example)
print(f"Is para-terphenyl: {result}, Reason: {reason}")