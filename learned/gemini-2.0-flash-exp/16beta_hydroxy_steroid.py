"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at the 16th carbon of a steroid core,
    where this hydroxy group has a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core SMARTS pattern (generic with implicit hydrogens)
    steroid_core_smarts = "[C]1[C][C]2[C]([C]1)([C][C]3[C]4[C]([C]2)[C]([C]3)([C])[C][C]4)[C]"

    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # Find substructure matches for steroid core
    core_match = mol.GetSubstructMatch(steroid_core_pattern)

    if not core_match:
        return False, "No steroid core found"

    # Map the match to the steroid core
    core_atoms = {}
    for i, atom_idx in enumerate(core_match):
       core_atoms[i+1] = atom_idx  # 1-based indexing for core atoms
    
    # Steroid ring carbon numbers with respect to atoms in core pattern match:
    # carbon 1: index 0
    # carbon 2: index 1
    # carbon 3: index 2
    # carbon 4: index 10
    # carbon 5: index 9
    # carbon 6: index 8
    # carbon 7: index 7
    # carbon 8: index 6
    # carbon 9: index 3
    # carbon 10: index 11
    # carbon 11: index 12
    # carbon 12: index 13
    # carbon 13: index 4
    # carbon 14: index 5
    # carbon 15: index 14
    # carbon 16: index 15

    # Get the atom index for the C16 carbon
    c16_idx = core_atoms.get(16)
    if c16_idx is None:
         return False, "Cannot determine index of C16"

    # Get the C16 atom object
    c16_atom = mol.GetAtomWithIdx(c16_idx)

    # Check for a hydroxyl group directly attached to C16
    has_oh_group = False
    for neighbor in c16_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen atom
            o_atom = neighbor
            # Check if O has one hydrogen attached
            h_count = 0
            for h_neighbor in o_atom.GetNeighbors():
                 if h_neighbor.GetAtomicNum() == 1: # Hydrogen
                     h_count += 1
            if h_count==1:
                has_oh_group = True
                o_c16_bond = mol.GetBondBetweenAtoms(c16_atom.GetIdx(), o_atom.GetIdx())
                # Verify beta configuration, checking for / bond only.
                if o_c16_bond.GetStereo() == Chem.BondStereo.STEREOUP:  # beta-OH should be an upward bond
                    return True, "16beta-hydroxy steroid found"

    if not has_oh_group:
        return False, "No hydroxyl group found at carbon 16"

    
    return False, "Hydroxyl group at C16, but not in beta configuration"