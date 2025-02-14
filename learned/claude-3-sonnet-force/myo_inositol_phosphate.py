"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:27767 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol backbone pattern
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No myo-inositol backbone found"

    # Exclude glycerophospholipids and glycosphingolipids
    glycerol_pattern = Chem.MolFromSmarts("C(C)(C)(OC)")
    sphingosine_pattern = Chem.MolFromSmarts("C(N)(C)(C)")
    if mol.HasSubstructMatch(glycerol_pattern) or mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Molecule is a glycerophospholipid or glycosphingolipid"

    # Look for phosphate groups (-O-P(=O)(-O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("[OX1,OX2]P(=O)([OX1,OX2])[OX1,OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"

    # Count number of phosphate groups
    n_phosphates = len(phosphate_matches)

    # Check for at least one phosphate group attached to the inositol ring
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)OP(=O)([OX1,OX2])[OX1,OX2]")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No phosphate group attached to the inositol ring"

    # Optionally, you can add additional checks for chemical properties or molecular weight

    return True, f"Contains myo-inositol backbone with {n_phosphates} phosphate group(s)"