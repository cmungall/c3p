"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is a sugar alcohol with a linear polyol structure derivable from the reduction of an aldose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alditol structure: [CH2OH][CHOH]*[CH2OH]
    alditol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CH2X4]")

    # Check if the molecule matches the alditol pattern
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "No continuous chain matching alditol structure found"
    
    # Count the number of hydroxyl (OH) groups
    num_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Ensure the count of OH groups is consistent with an alditol
    if num_oh_groups < 3:
        return False, f"Insufficient number of hydroxyl groups: {num_oh_groups} found"

    # Optionally, check if the molecule forms a continuous acyclic polyol
    is_cyclic = Chem.rdMolOps.FindMolChiralCenters(mol, includeUnassigned=True)
    if is_cyclic:
        return False, "Molecule contains a ring structure, not consistent with an acyclic alditol"

    return True, "Matches structure of an alditol with necessary OH groups"