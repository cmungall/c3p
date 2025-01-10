"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    Galactosylceramides have a beta-D-galactosyl moiety attached to a ceramide backbone,
    which includes amide linkages and long fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for galactosylceramide
    # Refined beta-D-galactosyl pattern considering possible variations
    galactose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](CO)C(O)C1")
    # Amide linkage: accounts for typical ceramide amide structure with carbonyl group
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    # Expanded hydrocarbon chain representation using a range
    long_chain_pattern = Chem.MolFromSmarts("C{8,}")

    # Match galactose
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "Missing beta-D-galactosyl group"
    
    # Match amide group (part of the ceramide backbone)
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage characteristic of ceramide"

    # Check for long-chain features indicative of fatty acid tails and sphingosine
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long-chain characteristic of ceramide"

    return True, "Contains structural motifs of a galactosylceramide: beta-D-galactosyl and ceramide backbone"