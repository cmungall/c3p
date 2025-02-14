"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2,
    consisting entirely of hydrogen atoms and saturated carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for atoms other than carbon and hydrogen
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 6 and atomic_num != 1:
            return False, f"Contains atom other than C and H: {atom.GetSymbol()}"

    # Check for rings (molecule must be acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, molecule is cyclic"

    # Check all bonds are single bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, "Contains multiple bonds, not all bonds are single"

    # Check all carbons are saturated (total valence of 4)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetTotalValence() != 4:
                return False, "Carbon atom not saturated (total valence not equal to 4)"

    # Calculate number of carbons and hydrogens
    num_carbons = 0
    num_hydrogens = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            num_carbons += 1
            num_hydrogens += atom.GetTotalNumHs()
        elif atom.GetAtomicNum() == 1:
            num_hydrogens += 1  # Count explicit hydrogens if any

    # Check general formula CnH2n+2
    if num_hydrogens != (2 * num_carbons + 2):
        return False, f"Does not fit general formula CnH2n+2 (found C{num_carbons}H{num_hydrogens})"

    return True, "Molecule is an acyclic saturated hydrocarbon fitting CnH2n+2"

__metadata__ = {   'chemical_class': {   'name': 'alkane',
                              'definition': 'An acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2, consisting entirely of hydrogen atoms and saturated carbon atoms.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet'},
        'message': None,
        'attempt': 1,
        'success': True}