"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    A ubiquinone is defined as a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core ubiquinone pattern (2,3-dimethoxy-5-methylbenzoquinone)
    core_smarts = "COc1c(C)c(=O)c(CO)c(=O)c1OC"  # Simplified SMARTS for core structure
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Unable to create core pattern"

    # Check for the core pattern in the molecule
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 2,3-dimethoxy-5-methylbenzoquinone not found"

    # Find the matching atoms for the core pattern
    matches = mol.GetSubstructMatches(core_pattern)
    if not matches:
        return False, "Core pattern not matched"

    # Assume position 6 is the carbon attached to the ring (not part of methoxy or methyl groups)
    # Identify attachment point for side chain (position 6)
    # In the core pattern, identify the atom that can have a side chain

    # For simplicity, define SMARTS pattern with labeled atom at position 6
    core_smarts_labeled = "COc1c(C)c(=O)c([C:6])c(=O)c1OC"
    core_pattern_labeled = Chem.MolFromSmarts(core_smarts_labeled)
    if core_pattern_labeled is None:
        return False, "Unable to create labeled core pattern"

    # Find matches with labeled atoms
    match_list = mol.GetSubstructMatch(core_pattern_labeled)
    if not match_list:
        return False, "Labeled core pattern not matched"

    # Get the atom index for position 6
    match_map = core_pattern_labeled.GetSubstructMatch(core_pattern_labeled)
    position_6_idx = None
    for atom in core_pattern_labeled.GetAtoms():
        if atom.HasProp('molAtomMapNumber') and atom.GetProp('molAtomMapNumber') == '6':
            core_atom_idx = atom.GetIdx()
            position_6_idx = match_list[core_atom_idx]
            break

    if position_6_idx is None:
        return False, "Could not locate position 6 in the molecule"

    # Get the atom at position 6
    position_6_atom = mol.GetAtomWithIdx(position_6_idx)

    # Check if there is a side chain attached at position 6
    side_chain_bonds = []
    for bond in position_6_atom.GetBonds():
        neighbor = bond.GetOtherAtom(position_6_atom)
        if neighbor.GetIdx() not in match_list:
            side_chain_bonds.append(bond)

    if not side_chain_bonds:
        return False, "No side chain attached at position 6"

    # Extract the side chain atoms
    side_chain_atoms = set()
    atoms_to_check = [position_6_atom]
    atoms_checked = set(match_list)
    while atoms_to_check:
        current_atom = atoms_to_check.pop()
        atoms_checked.add(current_atom.GetIdx())
        for bond in current_atom.GetBonds():
            neighbor = bond.GetOtherAtom(current_atom)
            if neighbor.GetIdx() not in atoms_checked:
                side_chain_atoms.add(neighbor.GetIdx())
                atoms_to_check.append(neighbor)

    if not side_chain_atoms:
        return False, "Side chain at position 6 is empty"

    # Create side chain molecule
    side_chain_mol = Chem.PathToSubmol(mol, list(side_chain_atoms))
    if side_chain_mol is None:
        return False, "Unable to isolate side chain"

    # Analyze the side chain
    # Check for sufficient length (number of carbons)
    num_carbons = sum(1 for atom in side_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 5:
        return False, "Side chain too short to be polyprenoid"

    # Check for presence of double bonds
    num_double_bonds = sum(1 for bond in side_chain_mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 1:
        return False, "Side chain lacks double bonds characteristic of polyprenoids"

    # Optionally, check for isoprene units (C5H8 units) by identifying patterns
    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
    isoprene_matches = side_chain_mol.GetSubstructMatches(isoprene_unit)
    if len(isoprene_matches) < 1:
        return False, "Side chain does not contain isoprene units"

    return True, "Contains 2,3-dimethoxy-5-methylbenzoquinone core with polyprenoid side chain at position 6"

__metadata__ = {   
    'chemical_class': {   
        'name': 'ubiquinones',
        'definition': 'Any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone; one of a group of naturally occurring homologues. The redox-active quinoid moiety usually carries a polyprenoid side chain at position 6, the number of isoprenoid units in which is species-specific. Ubiquinones are involved in the control of mitochondrial electron transport, and are also potent anti-oxidants.',
    }
}