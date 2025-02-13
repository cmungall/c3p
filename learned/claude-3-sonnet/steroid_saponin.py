"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:50047 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone (4 fused rings, 3 rings with 6 carbons, 1 ring with 5 carbons)
    steroid_pattern = Chem.MolFromSmarts("[C&r5,r6]1[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6]1[C&r5,r6]2[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r6]2[C&r5,r6]3[C&r5,r6][C&r5,r6][C&r5,r6][C&r5]3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for OH groups (hydroxysteroid)
    n_oh = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() < 2 for bond in atom.GetBonds()) > 0)
    if n_oh < 1:
        return False, "No hydroxyl groups found (not a hydroxysteroid)"
    
    # Check for glycosidic bonds (saponin)
    glycosidic_pattern = Chem.MolFromSmarts("[O-;!$(*=,#[!$(*=-,#[!$(*=-,-,-)])>1*]=C")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bonds found (not a saponin)"
    
    # Typical molecular weight range for steroid saponins
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 2000:
        return False, f"Molecular weight {mol_wt:.0f} Da outside typical range for steroid saponins"
    
    return True, "Contains steroid backbone with hydroxyl groups and glycosidic bonds"