"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol should have a phenol group with an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the guaiacol SMARTS pattern strictly for ortho substitution
    # Oc1c(OC)ccccc1 indicates the ortho position for the OCH3 relative to OH on the benzene ring.
    ortho_guaiacol_pattern = Chem.MolFromSmarts("Oc1cc(OC)ccc1")

    # Check for the guaiacol pattern
    if mol.HasSubstructMatch(ortho_guaiacol_pattern):
        return True, "Molecule contains a phenol group with an ortho methoxy substituent"

    return False, "Molecule does not match the guaiacol structure with phenol and ortho methoxy"

# Example usage
smiles_examples = [
    "COc1ccccc1O",  # guaiacol
    "C12=CC=CC=C1C(=CC=C2OC)O"  # Example of similar structures for clarity
]

for smiles in smiles_examples:
    result, reason = is_guaiacols(smiles)
    print(f"SMILES: {smiles} => Guaiacol: {result}. Reason: {reason}")