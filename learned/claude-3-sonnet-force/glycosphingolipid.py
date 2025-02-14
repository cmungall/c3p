"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: CHEBI:18095 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
    It is understood that the carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sphingoid/ceramide backbone
    ceramide_pattern = Chem.MolFromSmarts("[N;X3;H2,H1;!$(N(*)-*=[N,O,S])]")
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No sphingoid/ceramide backbone found"
    
    # Look for glycosidic linkage to the O-1 of the sphingoid
    for atom_idx in ceramide_matches:
        sphingoid_oxygen = mol.GetAtomWithIdx(atom_idx + 1)
        if sphingoid_oxygen.GetTotalNumHs() != 1:
            continue
        neighbors = mol.GetAtomWithIdx(atom_idx + 1).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
                glycosidic_oxygen = neighbor
                break
        else:
            continue
        glycosidic_neighbors = glycosidic_oxygen.GetNeighbors()
        for neighbor in glycosidic_neighbors:
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                return True, "Contains sphingoid/ceramide backbone with glycosidic linkage to O-1"

    return False, "No glycosidic linkage to O-1 of sphingoid/ceramide found"