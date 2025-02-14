"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: CHEBI:85556 2-monoglyceride

A monoglyceride in which the acyl substituent is located at position 2.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X3][CHX4][CH2X3]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester group (-O-C(=O)-) attached at position 2
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Find the carbon atom attached to the ester group
    ester_carbon = mol.GetAtomWithIdx(ester_matches[0][1]).GetNeighbors()[0].GetIdx()
    if ester_carbon != 2:  # Atom indices start from 0, so position 2 is index 2
        return False, "Ester group not attached at position 2 of glycerol"

    # Look for acyl chain (long carbon chain attached to ester)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_chain_matches) == 0:
        return False, "No acyl chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Acyl chain too short"

    return True, "Molecule has a glycerol backbone with an acyl chain attached at position 2"