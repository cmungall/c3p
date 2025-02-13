"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:36973 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdqueries

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound where monosaccharide units are joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycosidic bonds
    glycosidic_bond = Chem.MolFromSmarts("[OX2][CR]')
    if not mol.HasSubstructMatch(glycosidic_bond):
        return False, "No glycosidic bonds found"
    
    # Check for monosaccharide units
    monosaccharide_pattern = Chem.MolFromSmarts("[OX2][CR][CR][CR][CR][CR][CR]')
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    if not monosaccharide_matches:
        return False, "No monosaccharide units found"
    
    # Check for multiple monosaccharide units
    if len(monosaccharide_matches) < 2:
        return False, "Only one monosaccharide unit found, need at least two"
    
    # Check molecular weight range (typically 300-6000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 6000:
        return False, "Molecular weight outside typical range for oligosaccharides"
    
    # Count oxygen and carbon atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if o_count < c_count / 2:
        return False, "Too few oxygen atoms for an oligosaccharide"
    
    return True, "Contains multiple monosaccharide units joined by glycosidic linkages"