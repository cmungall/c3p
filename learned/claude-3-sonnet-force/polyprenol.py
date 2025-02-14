"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:27480 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has the general formula H-[CH2C(Me)=CHCH2]nOH, where the carbon skeleton
    is composed of more than one isoprene unit.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoprene unit pattern (C=C-C-C=C)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH2][CH2][CH2]=[CH]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Check for at least two isoprene units
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"
    
    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_match = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_match) != 1:
        return False, "Found more than one or no hydroxyl group"
    
    # Check for a hydrocarbon skeleton
    heteroatom_pattern = Chem.MolFromSmarts("[!#6;!#1]")  # Non-carbon, non-hydrogen atoms
    heteroatom_matches = mol.GetSubstructMatches(heteroatom_pattern)
    if len(heteroatom_matches) > 1:  # Excluding the terminal hydroxyl group
        return False, "Found heteroatoms other than the terminal hydroxyl group"
    
    return True, "Molecule has the general formula H-[CH2C(Me)=CHCH2]nOH with a hydrocarbon skeleton"

# Example usage
smiles = "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(=CCCC(=CCCC(=CCCC(=CCO)C)C)C)C)C)C)C)C"
print(is_polyprenol(smiles))