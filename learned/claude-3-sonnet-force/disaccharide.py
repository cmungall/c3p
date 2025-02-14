"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36973 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple

def is_disaccharide(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound where two monosaccharides are joined by a glycosidic bond.

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
    
    # Count number of monosaccharide rings
    sugar_rings = AllChem.MMRegionRings(mol)
    n_monosaccharides = len([r for r in sugar_rings if len(r) in [5, 6]])
    
    if n_monosaccharides != 2:
        return False, f"Found {n_monosaccharides} monosaccharide rings, expected 2"
    
    # Check for glycosidic bond (-O-)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond (-O-) found"
    
    # Check for other common disaccharide features
    has_pyranose_rings = any(len(r) == 6 for r in sugar_rings)
    has_furanose_rings = any(len(r) == 5 for r in sugar_rings)
    
    if has_pyranose_rings and has_furanose_rings:
        return True, "Contains two monosaccharide rings (pyranose and furanose) connected by a glycosidic bond"
    elif has_pyranose_rings:
        return True, "Contains two pyranose monosaccharide rings connected by a glycosidic bond"
    elif has_furanose_rings:
        return True, "Contains two furanose monosaccharide rings connected by a glycosidic bond"
    else:
        return False, "Unexpected ring structure for disaccharide"