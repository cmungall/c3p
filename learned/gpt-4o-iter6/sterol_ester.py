"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of 
    the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Sterol substructure pattern (modified for any typical sterol backbone)
    # using a more general cyclopentanoperhydrophenanthrene pattern
    sterol_pattern = Chem.MolFromSmarts("C1CCC2C3CC(C4CCC(C4)C3)C2C1")  # General sterol-like backbone
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Ester functional group pattern, general for esters
    # O-C(=O)- linkage indicative of ester bonds
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    return True, "Contains sterol backbone with ester linkage"

# Test the function with an example
example_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)OC(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
result, reason = is_sterol_ester(example_smiles)
print(f"Result: {result}, Reason: {reason}")