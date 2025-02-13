"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:36845 polyprenol

A polyprenol is defined as any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of more than one isoprene units.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.

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
    
    # Look for the isoprene unit pattern: CH2=C(C)CH=CH2
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C@H]([CH3])[CH2]=[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Polyprenols must have at least 2 isoprene units
    if len(isoprene_matches) < 2:
        return False, "Fewer than 2 isoprene units found"
    
    # Look for the terminal -OH group
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No terminal -OH group found"
    
    # Check for linear carbon chain
    linear_pattern = Chem.MolFromSmarts("[CH2X4]~[CH2X4]~[CH2X4]")
    linear_matches = mol.GetSubstructMatches(linear_pattern)
    if not linear_matches:
        return False, "Carbon chain is not linear"
    
    # Check for branching methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3X4]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) != len(isoprene_matches):
        return False, "Incorrect number of branching methyl groups"
    
    return True, "Molecule matches the polyprenol definition"