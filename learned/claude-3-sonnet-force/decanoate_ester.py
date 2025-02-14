"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:34892 decanoate ester
A fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for decanoate group (-C(=O)CCCCCCCCC)
    decanoate_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCCC")
    decanoate_matches = mol.GetSubstructMatches(decanoate_pattern)
    if len(decanoate_matches) == 0:
        return False, "No decanoate group found"

    # Look for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester linkage found"

    # Check if the decanoate group is part of the ester linkage
    for decanoate_idx in decanoate_matches:
        for ester_idx in ester_matches:
            if decanoate_idx[1] == ester_idx[0]:
                # Decanoate group is part of the ester linkage
                return True, "Contains a decanoate group linked via an ester bond"

    return False, "Decanoate group not linked via an ester bond"