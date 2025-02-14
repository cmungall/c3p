"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin has a carbon atom connected to both a hydroxyl group and a cyano group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a cyanohydrin core structure: C connected to -OH and -C#N
    cyanohydrin_core = Chem.MolFromSmarts("[CX4]([OH1])[C#N]")
    
    if mol.HasSubstructMatch(cyanohydrin_core):
      matches = mol.GetSubstructMatches(cyanohydrin_core)
      
      if len(matches) > 0 :
        return True, "Contains the core cyanohydrin substructure (carbon with OH and CN)"

    return False, "Does not contain the core cyanohydrin substructure"