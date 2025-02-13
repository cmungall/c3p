"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:51800 withanolide

A withanolide is any steroid lactone that is a C28 steroid with a modified side chain
forming a lactone ring and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a withanolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C28 steroid skeleton
    c28_pattern = Chem.MolFromSmarts("[C@]1(C)(CC[C@]2([H])[C@]3([H])CC[C@]4([H])C[C@@]([H])(O)[C@@H](O)C=C5[C@@]6([H])C[C@@H]7[C@]8([H])C[C@@H](O)[C@H](O)[C@]([H])(C9=CC(=O)O[C@H]9[C@@]1([H])C)[C@]7(C[C@H]8[C@]%10%11[C@@H]%12[C@@H]([C@H]([C@@H]%13[C@H]([C@@H]%14C)OC)C)OC%15)C)[C@@H]%16%17[C@H](C)[C@@H](O)[C@@H](C)[C@H](OC(=O)CC)[C@@H](C)[C@@H](OC)[C@@H](C)[C@@H](OC(=O)CC)[C@H](C)[C@H]3[C@]4([H])C[C@H]%18C%19C%20=C[C@H]([C@@H]6O5)[C@@H]2C")

    if not mol.HasSubstructMatch(c28_pattern):
        return False, "Does not contain C28 steroid skeleton"

    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts("[C&r5,r6]1[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6]1(=O)")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check molecular weight (typically 400-600 Da for withanolides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 600:
        return False, "Molecular weight outside typical range for withanolides"

    # Count number of oxygen atoms (typically 6-8 for withanolides)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6 or o_count > 8:
        return False, "Number of oxygen atoms outside typical range for withanolides"

    return True, "Contains C28 steroid skeleton with lactone ring and modified side chain"