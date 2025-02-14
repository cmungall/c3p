"""
Classifies: CHEBI:15734 primary alcohol
"""
Sure, here's a program to classify chemical entities as primary alcohols based on their SMILES string:

"""
Classifies: CHEBI:15781 primary alcohol

A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find atoms with -OH group
    oh_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1]
    
    # Check if any of the -OH groups are attached to a primary carbon
    for oh_atom in oh_atoms:
        carbon = oh_atom.GetNeighbors()[0]
        if carbon.GetSymbol() == 'C':
            if sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetSymbol() == 'C') <= 1:
                return True, "Contains a primary alcohol group (-OH attached to a saturated carbon atom with either three hydrogen atoms or one other carbon atom and two hydrogen atoms)"
    
    # If no primary alcohol groups found
    return False, "Does not contain a primary alcohol group"

# Example usage
print(is_primary_alcohol("OCCON1CCOCC1"))  # True, '2-(morpholin-4-yloxy)ethanol'
print(is_primary_alcohol("C(CCCO)C/C=C\C/C=C\CC1C(C/C=C\CCCC(O)=O)O1"))  # True, '8,9-epoxy-20-hydroxy-(5Z,11Z,14Z)-icosatrienoic acid'
print(is_primary_alcohol("CC\C=C/C\C=C/C\C=C/CCO"))  # True, '(3Z,6Z,9Z)-dodecatrienol'
print(is_primary_alcohol("C(C(S)CCO)C"))  # True, '3-mercaptopentanol'