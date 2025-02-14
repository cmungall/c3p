"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is not acyclic"
    
    # Count number of Hydroxyl groups ( -OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
         return False, "Molecule has fewer than 2 Hydroxyl groups"

    # Check for the presence of a carbon chain with hydroxyl groups
    polyol_pattern = Chem.MolFromSmarts("[CH2X4](O)-[CHX4](O)")
    polyol_matches = mol.GetSubstructMatches(polyol_pattern)
    if not polyol_matches:
        return False, "Molecule does not have a carbon chain with hydroxyl groups"

    # Check for allowed elements C, H, O only
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Molecule contains atoms other than C, H, and O"
    

    return True, "Molecule is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH"