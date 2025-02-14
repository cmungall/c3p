"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:36268 N-acylsphingosine

N-acylsphingosines are the parent compounds of the ceramide family, composed of
sphingosine having an unspecified fatty acyl group attached to the nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingosine backbone
    # 1. Double bond in the middle of a long carbon chain
    # 2. Primary alcohol (-CH2-OH) at one end
    # 3. Secondary amine (-NH-) at the other end
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H]")

    if (not mol.HasSubstructMatch(double_bond_pattern) or
        not mol.HasSubstructMatch(primary_alcohol_pattern) or
        not mol.HasSubstructMatch(secondary_amine_pattern)):
        return False, "No sphingosine backbone found"

    # Look for amide group (-N-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # Look for fatty acid chain (long carbon chain attached to amide)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern, maxMatches=1)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acid chain too short"

    return True, "Contains sphingosine backbone with fatty acyl chain attached via amide bond"