"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:38155 dodecanoate ester
A fatty acid ester in which the carboxylic acid component is lauric acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for lauric acid substructure (CCCCCCCCCCCCC(=O)O-)
    lauric_acid = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)[O-]")
    lauric_matches = mol.GetSubstructMatches(lauric_acid)
    if not lauric_matches:
        return False, "No lauric acid substructure found"

    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Check if ester is connected to lauric acid
    for ester_match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(ester_match[1])
        for nbr in ester_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 3 and nbr.GetTotalNumHs() == 0:
                if list(mol.GetAtomBonds(nbr.GetIdx())[2].GetBeginAtomIdx()) == lauric_matches[0]:
                    return True, "Contains lauric acid ester group"

    return False, "Lauric acid not connected to ester group"