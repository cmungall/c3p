"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is any aromatic amide obtained by acylation of aniline,
    meaning it has an amide group where the nitrogen is directly attached
    to a benzene ring (aromatic ring of six carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for anilide
    # Pattern explanation:
    # - [NX3H]: Amide nitrogen (three-coordinate nitrogen with one hydrogen)
    # - [CX3](=O): Carbonyl carbon (double-bonded to oxygen)
    # - [$(c1ccccc1)]: Benzene ring (six aromatic carbons in a ring)
    anilide_pattern = Chem.MolFromSmarts("[#6;a]1:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]:1[NX3H1][CX3](=O)[#6]")
    
    if mol.HasSubstructMatch(anilide_pattern):
        return True, "Contains an anilide moiety (amide nitrogen attached to benzene ring)"
    else:
        return False, "Does not contain an anilide moiety"