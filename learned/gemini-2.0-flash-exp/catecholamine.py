"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a benzene ring with hydroxyl groups at positions 1 and 2 (a catechol),
    and a 2-aminoethyl group attached at position 4, with possible substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for catechol moiety (1,2-dihydroxybenzene)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cccc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety (1,2-dihydroxybenzene) found"
    
    # Check for 2-aminoethyl group (-CH2-CH2-NH2)
    aminoethyl_pattern = Chem.MolFromSmarts("CCN")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No 2-aminoethyl group found"
        
    # Check that the 2-aminoethyl group is connected to the ring. We use a SMARTS pattern for that
    attachment_pattern = Chem.MolFromSmarts("c-C-C-N")
    if not mol.HasSubstructMatch(attachment_pattern):
        return False, "2-aminoethyl group not directly attached to the benzene ring"

    return True, "Contains a catechol moiety and a 2-aminoethyl group attached to the ring"