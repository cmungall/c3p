"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: CHEBI:37563 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins are a group of eight water-soluble vitamins that play important roles in cell metabolism,
    comprising of vitamin B1, B2, B3, B5, B6, B7, B9, and B12.

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

    # Check for common substructures and functional groups in B vitamins
    thiazole_pattern = Chem.MolFromSmarts("c1csc[n+]1")  # B1 (thiamine)
    if mol.HasSubstructMatch(thiazole_pattern):
        return True, "Contains thiazole ring, characteristic of vitamin B1 (thiamine)"

    isoalloxazine_pattern = Chem.MolFromSmarts("c1nc2nc3ccc(cc3nc2n1)C")  # B2 (riboflavin)
    if mol.HasSubstructMatch(isoalloxazine_pattern):
        return True, "Contains isoalloxazine ring system, characteristic of vitamin B2 (riboflavin)"

    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")  # B3, B5, B6
    pyridine_matches = mol.GetSubstructMatches(pyridine_pattern)
    if pyridine_matches:
        # Check for additional criteria for B3, B5, and B6
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH3+]")):  # B3 (niacin)
            return True, "Contains pyridine ring and ammonium group, characteristic of vitamin B3 (niacin)"
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)CC(O)=O")):  # B5 (pantothenic acid)
            return True, "Contains pyridine ring and acetic acid group, characteristic of vitamin B5 (pantothenic acid)"
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C=C")):  # B6 (pyridoxine)
            return True, "Contains pyridine ring and enone group, characteristic of vitamin B6 (pyridoxine)"

    pterin_pattern = Chem.MolFromSmarts("c1nc2nc(nc2n1)N")  # B9 (folate)
    if mol.HasSubstructMatch(pterin_pattern):
        return True, "Contains pterin ring system, characteristic of vitamin B9 (folate)"

    corrin_pattern = Chem.MolFromSmarts("[Co]")  # B12 (cobalamin)
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains corrin ring with cobalt, characteristic of vitamin B12 (cobalamin)"

    # Additional checks for molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1500:
        return False, "Molecular weight outside typical range for B vitamins"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too few rotatable bonds for a B vitamin"

    # No match found
    return False, "Structure does not match any known B vitamin substructure or criteria"