"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:25386 organohalogen compound
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for a carbon-halogen bond.
    # This pattern looks for a halogen connected to any atom by any bond, 
    # then checks if the connected atom is a carbon
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]~*")

    # Check if the molecule contains the pattern
    if mol.HasSubstructMatch(halogen_pattern):
        for match in mol.GetSubstructMatches(halogen_pattern):
            halogen_atom_index = match[0]
            neighbor_atom_index = match[1]
            
            halogen_atom = mol.GetAtomWithIdx(halogen_atom_index)
            neighbor_atom = mol.GetAtomWithIdx(neighbor_atom_index)

            if neighbor_atom.GetAtomicNum() == 6: # 6 is the atomic number for carbon
                return True, "Contains at least one carbon-halogen bond"

    return False, "No carbon-halogen bond found"