"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27279 vitamin D
Any member of a group of fat-soluble hydroxy seco-steroids that exhibit biological activity against vitamin D deficiency. 
Vitamin D can be obtained from sun exposure, food and supplements and is biologically inactive and converted into the biologically active calcitriol via double hydroxylation in the body.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Checks for a cholesterol-like steroid backbone with a characteristic side chain and hydroxylation pattern.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cholesterol-like steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@@]12[C@@]([C@H]([C@@H]3[C@H]([C@@H](C[C@@H]4[C@@]3(CC[C@@]4(O)[C@@H](O)[C@H]3[C@@]4([C@]2(C[C@@H](CC1)C)(CC[C@]34C)C)C)C)C)C)CC=C)C"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No cholesterol-like steroid backbone found"

    # Look for characteristic vitamin D side chain
    vit_d_pattern = Chem.MolFromSmarts("[C@@]12[C@H](C[C@H](C=C)C1)C[C@@H](O)[C@]2(O)CCC=C"
    if not mol.HasSubstructMatch(vit_d_pattern):
        return False, "Missing characteristic vitamin D side chain"

    # Check for hydroxylation pattern
    hydroxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIsAromatic() == False)
    if hydroxy_count < 1:
        return False, "No hydroxyl groups present"
    elif hydroxy_count == 1:
        return True, "Contains cholesterol-like steroid backbone with characteristic vitamin D side chain and one hydroxyl group"
    elif hydroxy_count == 2:
        return True, "Contains cholesterol-like steroid backbone with characteristic vitamin D side chain and two hydroxyl groups (calcitriol)"
    else:
        return True, f"Contains cholesterol-like steroid backbone with characteristic vitamin D side chain and {hydroxy_count} hydroxyl groups"

    return False, "Unable to classify molecule"