"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is characterized by a double bond between the 2nd and 3rd
    positions in the acyl chain and a Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for Coenzyme A
    coenzyme_a_pattern = Chem.MolFromSmarts("NCCSC(=O)")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "No Coenzyme A moiety found"

    # Define the SMARTS pattern for a double bond at the 2,3-position
    double_bond_2_3_pattern = Chem.MolFromSmarts("C=CC(=O)SCC")
    match_2_3_double_bond = mol.GetSubstructMatch(double_bond_2_3_pattern)

    if not match_2_3_double_bond:
        return False, "No double bond between the 2nd and 3rd positions"

    return True, "Contains a 2,3-double bond and Coenzyme A moiety"