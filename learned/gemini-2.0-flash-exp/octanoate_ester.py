"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to match octanoyl connected to an ester bond (C(=O)O)
    # Pattern matches a chain of 7 carbons followed by a carbonyl carbon connected to an oxygen via a single bond and that oxygen connected to any other atom.
    octanoyl_ester_pattern = Chem.MolFromSmarts("CCCCCCC[CX3](=O)[OX2][#6]")
    octanoyl_ester_matches = mol.GetSubstructMatches(octanoyl_ester_pattern)

     # SMARTS pattern to match ionized octanoyl connected to an ester bond
    octanoyl_ester_pattern_ion = Chem.MolFromSmarts("CCCCCCC[CX3](=O)[O-][#6]")
    octanoyl_ester_matches_ion = mol.GetSubstructMatches(octanoyl_ester_pattern_ion)

    if len(octanoyl_ester_matches) > 0 or len(octanoyl_ester_matches_ion) > 0 :
          return True, "Contains at least one octanoyl group connected to an ester"
    
    return False, "No octanoyl ester group found"