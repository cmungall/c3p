"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Glucosylceramides are cerebrosides with a glucose head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glucose moiety (pyranose ring with specific stereochemistry)
    # Note: This pattern matches beta-D-glucose specifically
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([OH0,OH1])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check for sphingoid base features:
    # - Long carbon chain
    # - NH group
    # - Two OH groups (or one OH in case of sphingosine)
    sphingoid_base = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2][CH]([OH1])[CH]([NH1])[CH2]O")
    if not mol.HasSubstructMatch(sphingoid_base):
        return False, "No sphingoid base found"

    # Check for amide bond (R-N-C(=O)-R)
    amide_pattern = Chem.MolFromSmarts("[NH1][C](=O)[CH2]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Check for fatty acid chain attached to amide
    # (at least 12 carbons typically)
    fatty_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No long fatty acid chain found"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Typical glucosylceramide should have:
    # - At least 30 carbons (glucose + sphingoid base + fatty acid)
    # - Around 8-9 oxygens (glucose has 6, plus 2-3 from ceramide)
    # - Exactly 1 nitrogen (from amide bond)
    if c_count < 30:
        return False, "Insufficient carbon atoms for glucosylceramide"
    if o_count < 8:
        return False, "Insufficient oxygen atoms for glucosylceramide"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom"

    return True, "Contains beta-D-glucose linked to ceramide with appropriate structural features"