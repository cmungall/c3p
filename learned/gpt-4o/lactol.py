"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by the intramolecular addition
    of a hydroxyl group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized lactol pattern:
    # A cyclic system where an oxygen atom (part of the ring, O1) is bonded to a carbon with a hydroxyl group
    lactol_pattern = Chem.MolFromSmarts('[O][C@H1,O1]([OH])')

    # Match the lactol pattern
    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal (lactol structure)"

    return False, "Lacks structural features of a lactol: cyclic hemiacetal"

# Test on identified lactol examples and review outcomes
examples = [
    "C/1(\C[C@H]2O[C@@H](C1)C[C@]3(C([C@H](C[C@@H](C[C@H](CC(O[C@@H]([C@@H](C)O)C[C@@]4(C\C(\[C@@H]([C@](C(C=C2)(C)C)(O4)O)OC(/C=C/C=C/CCC)=O)=C/C(OC)=O)[H])=O)O)O3)O)(C)C)O)=C\C(OC)=O",  # bryostatin 2
    "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose
    "O1C2=C(C(C[C@@]1(C3=CC(=C(C=C3)O)O)O)=O)C(=CC(=C2)O)O",  # (2S)-2-hydroxyeriodictyol
]

for smiles in examples:
    result, reason = is_lactol(smiles)
    print(f"SMILES: {smiles}\nIs Lactol: {result}\nReason: {reason}\n")