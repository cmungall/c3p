"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:24604 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecular entity that is a nucleoside
    in which one or more of the sugar hydroxy groups has been converted into a mono- or
    poly-phosphate. This includes both nucleotides and non-nucleotide nucleoside phosphates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nucleobase pattern (fused rings with N atoms)
    nucleobase_pattern = Chem.MolFromSmarts("c1ncnc2n1cncn2")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Look for sugar ring (5-membered ring with O and 4 C atoms)
    sugar_pattern = Chem.MolFromSmarts("C1OC(C)(C)C1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring found"

    # Look for phosphate groups (-O-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Check for glycosidic bond between nucleobase and sugar
    glycosidic_pattern = Chem.MolFromSmarts("c1ncnc2n1cncn2OC")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond between nucleobase and sugar"

    return True, "Contains a nucleobase, sugar ring, and phosphate group(s)"