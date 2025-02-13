"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies chemical entities as prenols based on their SMILES string.
Prenols are defined as 'Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).'
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenol(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX2H1]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No terminal alcohol group found"

    # Check for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("[CH2X4][CX3](=[CX3][CH2X4])[CH3]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene units found"

    # Count isoprene units
    n_isoprene = len(isoprene_matches)

    # Check for linear carbon skeleton
    linear_pattern = Chem.MolFromSmarts("[CH2X4][CH2X4]")
    linear_matches = mol.GetSubstructMatches(linear_pattern)
    if len(linear_matches) != n_isoprene:
        return False, "Carbon skeleton is not linear"

    # Check for double bonds in correct positions
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != n_isoprene:
        return False, "Incorrect number of double bonds"

    # Check for methyl groups in correct positions
    methyl_pattern = Chem.MolFromSmarts("[CX3][CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) != n_isoprene:
        return False, "Incorrect number of methyl groups"

    return True, "Molecule matches the prenol structural pattern"