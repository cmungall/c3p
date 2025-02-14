"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is a polyphenolic compound with glucoside linkages and possible catechol or pyrogallol residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for multiple phenolic rings
    phenol_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) < 2:
        return False, "Fewer than two phenolic rings found"

    # Look for catechol or pyrogallol moieties
    catechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)c(O)c1")
    if (not mol.HasSubstructMatch(catechol_pattern) and 
        not mol.HasSubstructMatch(pyrogallol_pattern)):
        return False, "No catechol or pyrogallol moieties found"

    # Check for glucose-like connections (typically as ethers or esters)
    glucose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No glucose or sugar connections found"

    # Additional checks such as molecular size can be considered,
    # but are omitted here for simplicity

    return True, "Contains polyphenolic structure with catechol/pyrogallol units and glucoside-like connections"