"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzoquinone core (1,4 or 1,2)
    benzoquinone_14_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")
    benzoquinone_12_pattern = Chem.MolFromSmarts("C1(=O)C(=O)C=CC=C1")
    #Check for naphthoquinone core (1,4 or 1,2)
    naphthoquinone_14_pattern = Chem.MolFromSmarts("C1(=O)C=CC2=CC=CC=C2C1=O")
    naphthoquinone_12_pattern = Chem.MolFromSmarts("C1(=O)C2=CC=CC=C2C(=O)C=C1")
    
    if not (mol.HasSubstructMatch(benzoquinone_14_pattern) or mol.HasSubstructMatch(benzoquinone_12_pattern) or
            mol.HasSubstructMatch(naphthoquinone_14_pattern) or mol.HasSubstructMatch(naphthoquinone_12_pattern)):
        return False, "No quinone core found"
        
    # Check for isoprenyl chain (at least two isoprene units)
    isoprenyl_pattern = Chem.MolFromSmarts("CC=CC")  # Basic isoprene unit
    isoprenyl_matches = mol.GetSubstructMatches(isoprenyl_pattern)
    if len(isoprenyl_matches) < 2:
        return False, "Less than two isoprenyl units found"


    return True, "Contains a quinone core and a polyprenyl side chain"