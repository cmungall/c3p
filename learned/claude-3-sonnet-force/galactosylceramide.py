"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:27798 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a ceramide with a galactose monosaccharide head group.

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
    
    # Check for sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[#8]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6]-[#6]-[#6](-[#6](-[#6]-[#6])=[#8])-[#6])=[#6]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Check for galactose head group
    galactose_pattern = Chem.MolFromSmarts("[OX2]C[CX4]([CX4]([CX4]([CX4]([OX2]-[#6])-[OX2])-[OX2])-[OX2])-[OX2]")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"
    
    # Check for amide bond linking backbone and head group
    amide_pattern = Chem.MolFromSmarts("[NX3](=[OX1])([#6])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond linking backbone and head group"
    
    # Check for long alkyl chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chains"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:
        return False, "Too few carbons for galactosylceramide"
    if o_count < 8:
        return False, "Too few oxygens for galactosylceramide"
    
    return True, "Contains sphingosine backbone with galactose head group and fatty acid chains"