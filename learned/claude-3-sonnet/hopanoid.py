"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:36806 hopanoid
A triterpenoid based on a hopane skeleton.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

def is_hopanoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hopane backbone pattern
    hopane_pattern = Chem.MolFromSmarts("[C@]12[C@@]([C@@](CC1)(C)[C@]3([C@@]4(CCCC(C)(C)[C@]4([H])CC3)C)CC2)C")
    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains a hopane skeleton"
    else:
        return False, "No hopane skeleton found"


# Example usage
smiles_list = [
    "CC(C)(O)[C@H]1CC[C@@]2(C)[C@H]1CC[C@]1(CO)[C@@H]2C[C@@H](O)[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3[C@H](O)C[C@@]12C",
    "CC(=C)[C@H]1CC[C@@]2(C)[C@H]1CC[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3C[C@H](O)[C@@]12C",
    "O=C(O)[C@H]1O[C@H](OCC(O)C(O)C(O)CC[C@@H](C2C3[C@](C4[C@@]([C@]5(C([C@@]6(C(C(C[C@H](C6)C)(C)C)CC5)C)CC4)C)(C)CC3)(C)CC2)C)[C@@H](O)[C@@H]([C@@H]1O)O",
    "CC(COC(C)=O)[C@H]1CC[C@@]2(C)[C@H]1CC[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C",
    "O=C1OC2C(O)[C@@H](OC1[C@@H]2O)OCC(O)C(O)C(O)CC[C@@H](C3C4[C@](C5[C@@]([C@]6(C([C@@]7(C(C(CCC7)(C)C)CC6)C)CC5)C)(C)CC4)(C)CC3)C",
    "CC(=C)[C@H]1CC[C@@]2(C)[C@H]1C[C@H](O)[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3C[C@H](O)[C@@]12C",
    "C([C@H]([C@@H](O)[C@H](CO)O)O)C[C@H]([C@]1(CC[C@@]2([C@]3(CC[C@@]4([C@]5(CCCC([C@@]5(CC[C@]4([C@@]3(CC[C@@]12[H])C)C)[H])(C)C)C)[H])[H])C)[H])C",
    "OC([C@@H]1C2[C@](C3[C@@]([C@]4(C([C@@]5(C(C(C[C@H](C5)C)(C)C)CC4)C)CC3)C)(C)CC2)(C)CC1)(C)C",
    "O(C1C(O)(C(O)C(C1N)O)CO)CC(O)C(O)C(O)C(O)CC([C@@H]2C3[C@](C4C=CC5[C@@]6(C(C([C@@H](C)CC6)(C)C)CC[C@]5([C@]4(C)CC3)C)C)(C)CC2)C",
    "OC1C(C2C(C3C(C4(C(C5(C(CC4)C(CC5)C(C)=C)C)CC3)C)(CC2)C)(CC1)C)(C)C"
]

for smi in smiles_list:
    print(smi)
    is_hop, reason = is_hopanoid(smi)
    print(f"Is hopanoid? {is_hop} - {reason}")
    print()