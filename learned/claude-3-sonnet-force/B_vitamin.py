"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: CHEBI:33219 B vitamins
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins are a group of water-soluble vitamins that play important roles in cell metabolism,
    including vitamins B1, B2, B3, B5, B6, B7, B9, and B12.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for common B vitamin substructures
    patterns = [
        # Thiamine (B1)
        Chem.MolFromSmarts("[C;H3][c;H1][c;H0][n;+]([C;H2][c;H1][s;H0]([C;H2][C;H2][O;H1][P;X4]=[O;X1])([C;H3])[C;H3])[c;H1][n;H1][c;H1]([C;H3])[C;H3]=[N;X2+]"),
        # Riboflavin (B2)
        Chem.MolFromSmarts("[n;X3]1[c;H0][n;H0][c;H1]2[c;H1]([n;X2][c;H1]1[n;X3][c;H0]2[C;X4][C;X4]=[O;X1])[C;X4][C;X4]=[O;X1]"),
        # Niacin (B3)
        Chem.MolFromSmarts("c1[c;H1][c;H1][c;H1][n;X2][c;H1]1[C;X3](=[O;X1])[O;H1]"),
        # Pantothenic acid (B5)
        Chem.MolFromSmarts("[C;H3][C;H1]([C;H3])([C;H3])[C;H2]([C;X3](=[O;X1])[N;X3][C;X4][C;X4][C;X3](=[O;X1])[O;H1])[O;H1]"),
        # Pyridoxine (B6)
        Chem.MolFromSmarts("[C;H3][c;H1]1[c;H0][c;H0][c;H1]([C;X3](=[O;X1])[O;H1])[c;H0][c;H0][n;X2]1[C;X4][O;H1]"),
        # Biotin (B7)
        Chem.MolFromSmarts("c1[c;H1][c;H1][c;H1][c;H1][c;H1]1[C;X3](=[O;X1])[N;X3][C;X4]1[C;X4][S;X2][C;X4]([C;H3])([C;H1])[N;X3][C;X3]1=[O;X1]"),
        # Folate (B9)
        Chem.MolFromSmarts("[c;H1]1[c;H1][c;H1][c;H1]([C;X3](=[O;X1])[N;X3][C;X4]([C;H2][C;H2][C;X3](=[O;X1])[O;H1])[C;X3](=[O;X1])[O;H1])[c;H1][c;H1]1"),
        # Cobalamin (B12)
        Chem.MolFromSmarts("[Co-3]1234[N;X4]5[C;X4]6=[N;+]1[C;X4](=[C;X3][C;X4]1=[N;+]2[C;X4](=[C;X3][C;X4]2=[N;+]3[C;X4]([C;H3])([C;H3])[C;X4]([C;H3])([C;H3])[C;H3][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]3([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]4([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]5([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]1([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]2([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;H3][C;X3](=[O;X1])[N;X3])[C;X4]6=[C;X3][C;X4]5=[N;+]1[C;X4](=[C;X3][C;X4]5=[N;+]6[C;X4]([C;H3])([C;H3])[C;X4]([C;H3])([C;H3])[C;H3][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]4([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;X4]2([C;H3])[C;X4][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;H3][C;X3](=[O;X1])[N;X3])[C;X4]4([C;H3])([C;H3])[C;X4][C;H3][C;X3](=[O;X1])[N;X3][C;X4]([C;H3])[C;H3][C;X3](=[O;X1])[N;X3])[C;H3]([C;H3])[C;H3]"),
    ]
    
    # Check if any pattern matches
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule matches a known B vitamin substructure"
    
    # Additional checks based on molecular properties
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 1000:
        return False, "Molecular weight too high for B vitamin"
    
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > 5:
        return False, "Too many rings for B vitamin"
    
    # Default to False if no pattern matched and properties didn't rule it out
    return False, "Molecule does not match any known B vitamin substructure"