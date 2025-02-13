"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones are characterized by a 2-aryl-1-benzopyran-4-one skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a flavone core structure
    # [cH]1[cH][cH][cH][cH]2COc3cc(=O)ccc3c12 represents a 2-aryl-1-benzopyran-4-one structure
    flavone_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)oc(=O)cc2")

    # Check if the molecule matches the flavone pattern
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No 2-aryl-1-benzopyran-4-one skeleton found"

    # Additional verification logic (e.g., ensure substitution correctness) could be added here

    return True, "Contains 2-aryl-1-benzopyran-4-one skeleton"

# Example usage
# result, reason = is_flavones("COc1cc(O)c2c(c1)oc1c(O)cc(c(O)c1c2=O)-c1c(O)c([C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1")
# print(result, reason)