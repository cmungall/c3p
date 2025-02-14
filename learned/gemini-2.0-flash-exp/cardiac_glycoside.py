"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are characterized by a steroid core, a lactone ring, and sugar residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str), True if it is a cardiac glycoside, False otherwise with a reason.
                If parsing fails return (None, None)
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Indicate an issue with SMILES parsing


    # 1. Check for steroid core
    #   SMARTS for basic steroid structure. Includes a fused 4-ring system.
    steroid_pattern = Chem.MolFromSmarts("[C]1[C]2[C][C]3[C]1[C]([C]([C]2)C)[C][C]4[C]3[C]([C]([C]4)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found."


    # 2. Check for Lactone ring. Often a 5-membered ring, attached to C17 of steroid
    lactone_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found."

    # 3. Check for sugar residues. This can be tricky, as sugars have multiple variations.
    #     Looking for a common sugar ring as a proxy
    sugar_pattern = Chem.MolFromSmarts("[O]1[CH][CH][CH]([O])[CH]([O])[CH]1")
    if not mol.HasSubstructMatch(sugar_pattern):
          return False, "No sugar residue found."



    # 4. Check for at least 2 OH groups. Note, this is a simplified check, may cause false positives
    oh_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalValence() == 2 and len([x for x in atom.GetNeighbors() if x.GetAtomicNum() == 1]) == 1 :
            oh_count += 1
    if oh_count < 2:
       return False, "Too few OH groups. Expected at least two."


    # Additional check on MW, as cardiac glycosides are typically larger.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 :
        return False, "Molecular weight too low to be a cardiac glycoside."

    return True, "Meets the criteria for a cardiac glycoside (steroid, lactone, sugar)."