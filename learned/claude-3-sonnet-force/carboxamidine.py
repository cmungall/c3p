"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:36562 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is a compound having the structure RC(=NR)NR2, containing
    the -C(=NH)NH2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all amidine groups (-C(=N)-N) in the molecule
    amidine_pattern = Chem.MolFromSmarts("C(=N-N)")
    amidine_matches = mol.GetSubstructMatches(amidine_pattern)
    
    # If no amidine groups found, it's not a carboxamidine
    if not amidine_matches:
        return False, "No carboxamidine (-C(=N)-N) group found"
    
    # Check if any of the amidine groups match the full carboxamidine pattern
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N-N)N")
    carboxamidine_matches = mol.GetSubstructMatches(carboxamidine_pattern)
    
    if carboxamidine_matches:
        return True, "Contains the carboxamidine (-C(=N)-N) group"
    else:
        return False, "Amidine group found but not a full carboxamidine structure"