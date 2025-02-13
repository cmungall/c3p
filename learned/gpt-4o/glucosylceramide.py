"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a cerebroside with a glucose head group and a sphingosine linked to a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-D-glucosyl pattern
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucosyl unit found"
    
    # Look for sphingosine backbone with amide bond (N-C(=O)) and long chains
    sphingosine_pattern = Chem.MolFromSmarts("NC(=O)C[C@@H](O)COC")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine structure with amide linkage found"
    
    # Look for a long aliphatic chain
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic chain (indicative of fatty acid) found"
    
    return True, "Contains glucose unit, sphingosine backbone with amide linkage, and long fatty acid chain"