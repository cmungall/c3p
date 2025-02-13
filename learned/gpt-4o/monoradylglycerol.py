"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("O[CH](CO)CO")
    
    # Check for the glycerol backbone match
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) != 1:
        return False, "Glycerol backbone not found or multiple present"

    # Identify and count valid substituents
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    alkyl_pattern = Chem.MolFromSmarts("C[C,C][C,C]")
    alk1enyl_pattern = Chem.MolFromSmarts("C=C")

    acyl_match_count = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_match_count = len(mol.GetSubstructMatches(alkyl_pattern))
    alk1enyl_match_count = len(mol.GetSubstructMatches(alk1enyl_pattern))

    # We expect exactly one substituent type to match
    if acyl_match_count == 1 or alkyl_match_count == 1 or alk1enyl_match_count == 1:
        return True, "Contains a glycerol backbone with a single acyl, alkyl or alk-1-enyl substituent"
    
    return False, "Glycerol backbone does not have exactly one suitable substituent"