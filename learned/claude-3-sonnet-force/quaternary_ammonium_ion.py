"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:53054 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define quaternary N+ pattern
    quat_n = rdqueries.AtomIsPositivelyCharged() & rdqueries.AtomIsAliphaticNitrogen()
    
    # Check for quaternary N+ with 4 substituents
    quat_n_matches = mol.GetSubstructMatches(quat_n)
    
    if not quat_n_matches:
        return False, "No quaternary nitrogen found"
    
    for match in quat_n_matches:
        atom = mol.GetAtomWithIdx(match)
        if sum(1 for _ in atom.GetNeighbors()) != 4:
            return False, "Quaternary nitrogen does not have 4 substituents"
    
    return True, "Contains a positively charged quaternary nitrogen with 4 substituents"