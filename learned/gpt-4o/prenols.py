"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is an alcohol consisting of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isoprene pattern: [CH2-C(Me)=CH-CH2], representing `C-C(C)=C-C`
    isoprene_pattern = Chem.MolFromSmarts("[CH2][C](C)=C[CH2]") 
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    # Check if there's at least one isoprene unit
    if len(isoprene_matches) < 1:
        return False, "No isoprene units found"

    # Verify the presence of a terminal alcohol group (-OH)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
        return False, "Missing terminal OH group"

    # Determine the connectivity of the alcohol group, typically terminal
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and len(list(atom.GetBonds())) == 1:  # Oxygen with only one bond
            bond = atom.GetBonds()[0]
            connected_atom = bond.GetOtherAtom(atom)
            # Look for isoprene unit connection
            if connected_atom.GetAtomicNum() == 6 and connected_atom.GetDegree() >= 1:
                return True, "Contains isoprene units with terminal alcohol group"
    
    return True, "Contains isoprene units with terminal alcohol group"