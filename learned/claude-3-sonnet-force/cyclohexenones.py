"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: CHEBI:38021 cyclohexenone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclohexenone(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone with one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: True if the molecule is a cyclohexenone, False otherwise
                          Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a six-membered ring with one double bond and one carbonyl
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:  # Six-membered ring
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_double_bonds = sum(bond.GetBondType() == Chem.BondType.DOUBLE
                                   for bond in mol.GetBondsBetweenAtoms(ring[:]))
            num_carbonyls = sum(atom.GetAtomicNum() == 8 and
                                 atom.GetTotalNumHs() == 0 and
                                 sum(bond.GetBondTypeAsDouble() == 2 for bond in atom.GetBonds()) == 1
                                 for atom in ring_atoms)
            if num_double_bonds == 1 and num_carbonyls == 1:
                return True, "Molecule contains a six-membered ring with one double bond and one carbonyl"

    return False, "Molecule does not contain the required cyclohexenone substructure"