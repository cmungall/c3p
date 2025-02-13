"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: CHEBI:28016 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of multiple biological molecule fragments
    biological_fragments = [
        Chem.MolFromSmarts("[NX3+]"),  # Amino groups (peptides/proteins)
        Chem.MolFromSmarts("[OX2H1]"),  # Hydroxyl groups (carbohydrates)
        Chem.MolFromSmarts("[cX2]1[nX2]c[nX2]c1"),  # Nucleic acid bases
        Chem.MolFromSmarts("[SX2H1]"),  # Thiol groups (glutathione, cysteine, etc.)
        Chem.MolFromSmarts("[CX4H1,CX3]([CX4H2])([CX4H2])[CX4H3]"),  # Fatty acid chains
    ]
    fragment_counts = [len(mol.GetSubstructMatches(frag)) for frag in biological_fragments]
    
    # Require at least 2 fragments from different classes
    unique_fragments = len([count for count in fragment_counts if count > 0])
    if unique_fragments < 2:
        return False, "Less than 2 distinct biological molecule fragments found"
    
    # Check for covalent linkage between fragments
    linker_pattern = Chem.MolFromSmarts("[#6]~[#6]")
    linker_matches = mol.GetSubstructMatches(linker_pattern)
    if not linker_matches:
        return False, "No covalent linkage found between fragments"
    
    return True, "Contains at least 2 distinct biological molecule fragments covalently linked together"