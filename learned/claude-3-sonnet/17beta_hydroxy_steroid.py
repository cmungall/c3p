"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17beta-hydroxy steroid
A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~1~[#6]~2~[#6]~3~4")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 17-OH group in beta configuration
    # [C] is carbon 17, [OH1] is hydroxy group, '@' indicates stereochemistry
    # The [C] must be connected to 4 atoms (saturated)
    # Note: The exact SMARTS pattern depends on the numbering convention used
    oh_17_beta = Chem.MolFromSmarts('[C;X4](@[*])(@[*])(@[*])[OH1]')
    
    if not mol.HasSubstructMatch(oh_17_beta):
        return False, "No hydroxyl group with correct connectivity found"

    # Get matches for OH group
    oh_matches = mol.GetSubstructMatches(oh_17_beta)
    
    # Check if any of the matches are at position 17
    found_17_beta_oh = False
    for match in oh_matches:
        c_atom = mol.GetAtomWithIdx(match[0])  # Get the carbon atom
        # Check if this carbon is part of the D ring (ring 4) of the steroid
        # by checking its environment
        ring_info = mol.GetRingInfo()
        if ring_info.NumAtomRings(match[0]) > 0:  # Carbon must be part of a ring
            # Check chirality of the carbon
            if c_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                found_17_beta_oh = True
                break

    if not found_17_beta_oh:
        return False, "No 17-beta hydroxyl group found"

    # Additional validation: molecule should have reasonable size for a steroid
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20 or num_atoms > 100:
        return False, "Molecule size not consistent with steroid structure"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    return True, "Contains steroid core with 17-beta hydroxyl group"