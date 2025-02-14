"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:36512 acyl-CoA
Acyl-CoA is defined as a thioester that results from the formal condensation of the thiol group
of coenzyme A with the carboxy group of any carboxylic acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("C1(N2C3=C(N=C(N=C3N)N)N=C2)OC(COP(=O)(O)OP(=O)(O)OCC(C(O)C(=O)NCCC(=O)NCCS)O)C(O)C1O")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No coenzyme A moiety found"

    # Look for thioester (-C(=O)-S-) group attached to CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Check if thioester is connected to CoA
    for match in thioester_matches:
        thioester_atom = mol.GetAtomWithIdx(match)
        if thioester_atom.GetNeighbors()[1].IsInRing():
            return True, "Contains coenzyme A moiety connected to a thioester group"

    return False, "Thioester group not connected to coenzyme A moiety"