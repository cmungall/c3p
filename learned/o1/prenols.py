"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH,
    in which the carbon skeleton is composed of one or more isoprene units.

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

    # Check for alcohol group (-OH)
    OH_pattern = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(OH_pattern):
        return False, "No alcohol group (-OH) found"

    # Define isoprene unit pattern: CH2-C(Me)=CH-CH2
    isoprene_pattern = Chem.MolFromSmarts('[CH2]-[C]([CH3])=C-[CH2]')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) == 0:
        return False, "No isoprene units found"

    # Check that isoprene units are connected head-to-tail
    # Build a graph of isoprene units
    isoprene_atoms = set()
    for match in isoprene_matches:
        isoprene_atoms.update(match)

    # Get all atoms in the molecule
    total_atoms = set(range(mol.GetNumAtoms()))

    # Remove isoprene atoms from total atoms, the remaining should be minimal
    remaining_atoms = total_atoms - isoprene_atoms
    if len(remaining_atoms) > 5: # Allowing for OH group and possible variations
        return False, "Isoprene units are not connected head-to-tail or extra atoms detected"

    return True, "Contains alcohol group and one or more isoprene units connected head-to-tail"