"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for rings
    if not mol.GetRingInfo().IsRingAtom(0): # if the first atom is not in a ring, fail fast
        if not any(mol.GetRingInfo().IsRingAtom(atom.GetIdx()) for atom in mol.GetAtoms()):
            return False, "Molecule does not contain a ring"

    # 2. Identify hemiacetal carbon (C-OH and C-O-C in a ring). Look for -C(OH)-O- pattern in a ring
    hemiacetal_pattern = Chem.MolFromSmarts("[CX4]([OX2H1])[OX2]1[#6]~[#6]1")
    matches = mol.GetSubstructMatches(hemiacetal_pattern)

    if not matches:
        hemiacetal_pattern_2 = Chem.MolFromSmarts("[CX4]([OX2H1])[OX2][#6]1~[#6]~[#6]1") # allow for other ring patterns
        matches = mol.GetSubstructMatches(hemiacetal_pattern_2)

    if not matches:
        hemiacetal_pattern_3 = Chem.MolFromSmarts("[#6][CX4]([OX2H1])[OX2]1[#6]~[#6]1") # allow for carbon before hemiacetal
        matches = mol.GetSubstructMatches(hemiacetal_pattern_3)

    if not matches:
        hemiacetal_pattern_4 = Chem.MolFromSmarts("[#6][CX4]([OX2H1])[OX2][#6]1~[#6]~[#6]1") # allow for carbon before hemiacetal
        matches = mol.GetSubstructMatches(hemiacetal_pattern_4)


    if not matches:
        return False, "No hemiacetal group found in ring system"

    # Additional check: Ensure that this hemiacetal is not part of a larger acetal or ketal
    # This is not perfect, but should prevent some false positives with structures that contain two ether oxygens attached to the central carbon.
    central_carbon_smarts = Chem.MolFromSmarts("[CX4]([OX2H1])[OX2]")
    
    for match in matches:
        central_carbon = mol.GetAtomWithIdx(match[0])
        connected_oxygens = 0
        for neighbor in central_carbon.GetNeighbors():
          if neighbor.GetAtomicNum() == 8:
            connected_oxygens += 1
        if connected_oxygens > 2:
           return False, "Hemiacetal carbon is connected to more than 2 oxygens."
   
    return True, "Contains a cyclic hemiacetal structure"