"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    The structure must feature a hydroxyl group at position 17 in a beta configuration on a steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic SMARTS pattern for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C=CC2)CC4=C3C=CC=C4C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"

    # Identify 17beta-hydroxy group pattern on steroid
    # The pattern explores beta-attachment on the typical steroid carbon backbone
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]1(C)C[C@@H]([C@H](O)[H])CC[C@@H](C)C1")
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() == False:
            idx = atom.GetIdx()
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(8) == 1:  # Hydroxyl group
                # Check stereochemistry
                if mol.GetAtomWithIdx(idx).GetChiralTag() == Chem.rdchem.CHI_TETRAHEDRAL_CCW:
                    return True, "17beta-Hydroxy steroid identified"
    
    return False, "Does not have a 17beta-hydroxy configuration on steroid scaffold"