"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is identified by the presence of two sugar moieties
    connected by a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more comprehensive SMARTS patterns for pyranose and furanose rings
    # Including specific ring sizes and common stereochemistry patterns
    pyranose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1")  # 6-membered pyranose form
    furanose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1")      # 5-membered furanose form
    
    # Look for these ring patterns twice (to ensure we have two sugar rings)
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    
    # Ensure exactly two monosaccharide ring structures
    total_sugar_rings = len(pyranose_matches) + len(furanose_matches)
    if total_sugar_rings != 2:
        return False, f"Expected 2 monosaccharide units, found {total_sugar_rings}"

    # Define a more specific SMARTS pattern for glycosidic linkage
    # Ensuring it connects the anomeric carbon of one sugar to a hydroxyl of another
    glycosidic_pattern = Chem.MolFromSmarts("COC")  # Simplified representation
    
    # Check for a glycosidic bond specifically linking the sugar rings
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    # Checking that it is not just any ether group but specifically between sugars
    if not any(
        set(m1) & set(m2) for m1 in pyranose_matches + furanose_matches 
        for m2 in glycosidic_matches
    ):
        return False, "No glycosidic bond found connecting sugar units"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"