"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:24996 octadecadienoic acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is any straight-chain C18 polyunsaturated fatty acid
    having two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Ensure the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight-chain fatty acid"

    # Ensure the molecule is unbranched
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetDegree() > 3:
                return False, "Molecule is branched; not a straight-chain fatty acid"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Number of carbon atoms is {c_count}, expected 18"

    # Count the number of carbon-carbon double bonds
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                cc_double_bonds += 1
    if cc_double_bonds != 2:
        return False, f"Number of C=C double bonds is {cc_double_bonds}, expected 2"

    return True, "Molecule is a straight-chain C18 fatty acid with two C=C double bonds"