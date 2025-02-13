"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:28167 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is a flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-arylchromen-4-one skeleton
    # 1. Find the chromone core
    chromone_pattern = Chem.MolFromSmarts("O=C1C=C(Oc2ccccc2)C=CC1")
    chromone_matches = mol.GetSubstructMatches(chromone_pattern)
    if not chromone_matches:
        return False, "No chromone core found"
    
    # 2. Check for aromatic ring attached to the chromone
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    for match in chromone_matches:
        chromone_atom = mol.GetAtomWithIdx(match[2])
        for neighbor in chromone_atom.GetNeighbors():
            if neighbor.IsInRingSize(6):
                aromatic_ring = Chem.FindMolRing(mol, thirdSeedIdx=neighbor.GetIdx())
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in aromatic_ring):
                    return True, "Contains 2-arylchromen-4-one skeleton"
    
    return False, "No 2-arylchromen-4-one skeleton found"