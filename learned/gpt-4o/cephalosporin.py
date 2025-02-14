"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    A cephalosporin has a beta-lactam structure with a six-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define beta-lactam (azetidinone) ring SMARTS pattern with four atoms: O=C1N(C=O)C1
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N(C)C1=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Define six-membered dihydrothiazine ring surrounded by a S atom SMARTS pattern
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1=C(SC2)C(N2)=C1")
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found"

    return True, "Contains beta-lactam and dihydrothiazine rings typical of cephalosporins"

# Example usage
smiles_examples = [
    "C1C(=O)N(C2=CCS(=O)(=O)C2C1=C(O)O)C(=O)O",  # Cephalosporin C variant
    "C1C(=O)N(C2=CCS(=O)(=O)C2C1=C(O)Cl)C(=O)O",  # Variant with a Cl substitution
    "C1C(=O)N(C2=CCN(S)C2C1=C(O)O)C(=O)O",       # Not a cephalosporin (no thiazine ring)
]

for smiles in smiles_examples:
    is_cef, reason = is_cephalosporin(smiles)
    print(f"SMILES: {smiles} - Is Cephalosporin? {is_cef}. Reason: {reason}")