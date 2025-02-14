"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon having one carbon-carbon double bond and the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acyclicity (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic"

    # Check that all atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Contains atom other than carbon or hydrogen: {atom.GetSymbol()}"

    # Count number of carbon-carbon double bonds
    num_cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_cc_double_bonds += 1

    if num_cc_double_bonds != 1:
        return False, f"Contains {num_cc_double_bonds} carbon-carbon double bonds, requires exactly 1"

    # Compute molecular formula and check if it fits CnH2n
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    if h_count != 2 * c_count:
        return False, f"Molecular formula does not fit CnH2n: C{c_count}H{h_count}"

    return True, "Molecule is an acyclic hydrocarbon with exactly one carbon-carbon double bond and fits CnH2n formula"

__metadata__ = {
    'chemical_class': {
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.'
    }
}