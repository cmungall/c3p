"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as any aromatic amide obtained by acylation of aniline.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Anilide feature patterns: Phenyl group directly bonded to a nitrogen in an amide linkage
    anilide_patterns = [
        # Direct phenyl linkage to amide nitrogen
        Chem.MolFromSmarts("c1ccccc1N(C(=O)[#6])"),
        # Allowing for substitutions on the phenyl ring
        Chem.MolFromSmarts("c1ccc(cc1)N(C(=O)[#6])"),
        # Allowing for substitutions on the amine nitrogen
        Chem.MolFromSmarts("c1ccccc1N([#0,#6,#7,#8])C(=O)[#6]"),
        # More generalized pattern to cover diverse anilide configurations
        Chem.MolFromSmarts("c1ccc(cc1)N([#0,#6,#7,#8])C(=O)[#6]")
    ]
    
    # Check if any of the patterns match, indicating an anilide structure
    for pattern in anilide_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an aniline-like aromatic amide structure"

    return False, "Does not match anilide structural criteria"