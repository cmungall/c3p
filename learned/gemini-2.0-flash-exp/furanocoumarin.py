"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is a coumarin core with a fused furan ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for coumarin and furan
    coumarin_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)OC=C") #Basic coumarin
    furan_pattern = Chem.MolFromSmarts("c1occc1")  # Furan ring
    
    # Check for coumarin core
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin core found"

    # Check for furan ring
    if not mol.HasSubstructMatch(furan_pattern):
          return False, "No furan ring found"

    # Check for the fusion
    # The following check is not ideal, since it does not guarantee that both rings are connected
    # ideally we should be looking for the full fragment.
    # but it was too difficult to generate a smarts for this.
    # this check is only for a minimum overlap and should not cause too many false positives.
    
    combined_pattern = Chem.MolFromSmarts("[c1ccccc1C(=O)OC=C]~[c1occc1]")
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "Furan ring is not fused to the coumarin core"

    return True, "Contains a coumarin core with a fused furan ring"