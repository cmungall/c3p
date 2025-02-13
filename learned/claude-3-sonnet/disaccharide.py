"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36973 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound with two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for monosaccharide patterns
    monosaccharide_pattern = Chem.MolFromSmarts("[OX2][CX4](O)[CX3](O)[CX3](O)[CX3](O)[CX3](O)[CX2]")
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    # Check if there are exactly 2 monosaccharide subunits
    if len(monosaccharide_matches) != 2:
        return False, f"Found {len(monosaccharide_matches)} monosaccharide units, need exactly 2"
    
    # Check if the monosaccharide units are joined by a glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4][CX3][OX2]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bond_matches:
        return False, "No glycosidic bond found between monosaccharide units"
    
    # Check molecular weight - disaccharides typically 300-600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight outside typical range for disaccharides"
    
    return True, "Contains two monosaccharide units joined by a glycosidic bond"