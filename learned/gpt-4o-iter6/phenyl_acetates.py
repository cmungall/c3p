"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by condensation of acetic acid with
    the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for acetoxy group connecting to an aromatic ring
    acetoxy_aromatic_pattern = Chem.MolFromSmarts("cCOC(=O)C")
    if mol.HasSubstructMatch(acetoxy_aromatic_pattern):
        return True, "Molecule is a phenyl acetate with the acetate ester linked to an aromatic carbon"

    # As an additional coverage for the acetate bound on catechol/phenol derivatives
    acetoxy_phenol_pattern = Chem.MolFromSmarts("OcCOC(=O)C")
    if mol.HasSubstructMatch(acetoxy_phenol_pattern):
        return True, "Molecule features phenolic oxygen linked to phenyl acetate"

    return False, "Molecule does not exhibit key phenyl acetate structural traits"