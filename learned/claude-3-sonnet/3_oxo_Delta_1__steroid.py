"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    These are steroids with a ketone at position 3 and a double bond between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid ring system (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) that's part of the first ring of the steroid
    ketone_pattern = Chem.MolFromSmarts("[#6]~1~[#6]=,:[#6]~[#6](=[O])~[#6]~[#6]~[#6]~2")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone found at position 3"

    # Check for double bond between positions 1 and 2
    # This pattern specifically looks for the double bond between C1 and C2,
    # ensuring it's part of the ring with the ketone at C3
    double_bond_pattern = Chem.MolFromSmarts("[#6]~1~[#6]=,:[#6]~[#6](=[O])~[#6]~[#6]~[#6]~2")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond between positions 1 and 2"

    # Additional check to ensure the molecule has reasonable size for a steroid
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Steroids typically have more than 20 atoms
        return False, "Molecule too small to be a steroid"

    return True, "Molecule contains steroid core with 3-oxo group and Delta(1) double bond"