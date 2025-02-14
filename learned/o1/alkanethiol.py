"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thiol pattern: sulfur (S) with one hydrogen (H) attached and connected to carbon (C)
    thiol_pattern = Chem.MolFromSmarts("[#16H1]-[#6]")
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "Does not contain an -SH group attached to an alkyl group"

    # Define amide bond pattern: carbonyl group attached to nitrogen (peptide bond)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bonds (peptide), not an alkanethiol"

    return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"