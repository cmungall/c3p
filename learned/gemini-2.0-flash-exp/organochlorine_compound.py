"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for C-Cl bond
    for bond in mol.GetBonds():
      atom1 = bond.GetBeginAtom()
      atom2 = bond.GetEndAtom()
      if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 17) or \
         (atom1.GetAtomicNum() == 17 and atom2.GetAtomicNum() == 6):
        return True, "Contains at least one carbon-chlorine bond"
    
    return False, "Does not contain any carbon-chlorine bond"