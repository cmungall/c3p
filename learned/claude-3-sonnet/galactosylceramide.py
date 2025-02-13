"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide has a galactose sugar connected to a ceramide 
    (sphingosine + fatty acid) via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for galactose pattern (6-membered ring with 4 OH groups and 1 CH2OH)
    # Note: This is a simplified pattern that matches pyranose sugars
    sugar_pattern = Chem.MolFromSmarts("[CH2][CH]1[CH]([OH])[CH]([OH])[CH]([OH])[CH]([OH])O1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No galactose moiety found"

    # Look for amide bond (-NH-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NH][C](=[O])[CH2,CH]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Look for long carbon chains (sphingosine base and fatty acid)
    carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")
    chain_matches = mol.GetSubstructMatches(carbon_chain)
    if len(chain_matches) < 2:
        return False, "Missing required long carbon chains"

    # Look for glycosidic bond pattern (C-O-C connecting sugar to sphingosine)
    glycosidic_pattern = Chem.MolFromSmarts("[CH2]O[CH]1O[CH][CH][CH][CH][CH]1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"

    # Check for characteristic OH groups on sphingosine
    sphingosine_oh_pattern = Chem.MolFromSmarts("[CH]([OH])[CH]([NH])")
    if not mol.HasSubstructMatch(sphingosine_oh_pattern):
        return False, "Missing characteristic OH groups on sphingosine"

    # Count carbons and oxygens to ensure molecule is large enough
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 24:  # Typical galactosylceramides have at least 24 carbons
        return False, "Too few carbons for galactosylceramide"
    if o_count < 7:  # At least 7 oxygens (4 from galactose, 2 from sphingosine, 1 from amide)
        return False, "Too few oxygens for galactosylceramide"

    # Calculate molecular weight - should be substantial due to long chains
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for galactosylceramide"

    # Optional: check for sulfate groups which are present in some examples
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    base_reason = "Contains galactose connected to ceramide (sphingosine + fatty acid) via glycosidic bond"
    if has_sulfate:
        return True, base_reason + " with sulfate modification"
    return True, base_reason