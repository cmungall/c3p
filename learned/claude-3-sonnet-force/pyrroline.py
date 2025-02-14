"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:22990 pyrroline
Any organic heteromonocyclic compound with a structure based on a dihydropyrrole.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a structure based on a dihydropyrrole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a single ring with 5 atoms, 4 carbon and 1 nitrogen
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring found"
    if len(ring_info.AtomRings()[0]) != 5:
        return False, "Ring size is not 5"
    ring_atoms = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring_info.AtomRings()[0]]
    if ring_atoms.count(6) != 4 or ring_atoms.count(7) != 1:
        return False, "Ring does not contain 4 carbon and 1 nitrogen atoms"

    # Check for double bond in the ring
    bond_types = [mol.GetBondBetweenAtoms(ring_info.AtomRings()[0][i], ring_info.AtomRings()[0][(i+1)%5]).GetBondType() for i in range(5)]
    if Chem.rdchem.BondType.DOUBLE not in bond_types:
        return False, "No double bond found in the ring"

    # Check for exocyclic substituents
    allowed_substituents = [Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP2]
    for atom_idx in ring_info.AtomRings()[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_info.AtomRings()[0]:
                if neighbor.GetHybridization() not in allowed_substituents:
                    return False, "Substituent is not allowed for pyrrolines"

    return True, "Molecule contains a 5-membered ring with 4 carbon and 1 nitrogen atoms, a double bond, and allowed substituents"