"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule can be classified as a potential glycosaminoglycan 
    based on its SMILES string. Glycosaminoglycans are polysaccharides with a 
    substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a potential glycosaminoglycan, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("O[C;R]O")

    # Define pattern for aminomonosaccharide sugars with nitrogen in a carbon ring
    amino_sugar_pattern = Chem.MolFromSmarts("[C;R][NX3;H2,H1]")

    # Detect glycosidic-type linkages
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Detect amino sugar residues
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)
    
    # Determine classification based on counts
    if len(glycosidic_matches) > 0 and len(amino_sugar_matches) > 0:
        # Considering proportion of amino sugars to polysaccharide units
        if len(amino_sugar_matches) / max(1, len(glycosidic_matches)) > 0.5:  # Adjusted threshold
            return True, "Contains glycosidic structures with a substantial proportion of amino sugar units"
        else:
            return False, "Insufficient proportion of amino sugar units"
    
    if len(glycosidic_matches) == 0:
        return False, "Lacks characteristic glycosidic linkages for polysaccharides"
    if len(amino_sugar_matches) == 0:
        return False, "Lacks characteristic amino sugar residues"

    return None, "Could not be definitively classified as a glycosaminoglycan"