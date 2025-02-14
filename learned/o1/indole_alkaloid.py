"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24842 indole alkaloid
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms (since alkaloids contain nitrogen)
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if 7 not in atom_nums:
        return False, "No nitrogen atoms found (not an alkaloid)"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"

    # Store ring properties
    rings = []
    for ring_atoms in atom_rings:
        ring = {
            'atoms': ring_atoms,
            'size': len(ring_atoms),
            'aromatic': all([mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms]),
            'has_nitrogen': any([mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring_atoms])
        }
        rings.append(ring)

    # Search for fused rings
    found_indole = False
    num_rings = len(rings)
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # Check if rings are fused (share at least two atoms)
            set_i = set(rings[i]['atoms'])
            set_j = set(rings[j]['atoms'])
            shared_atoms = set_i.intersection(set_j)
            if len(shared_atoms) >= 2:
                # Check for indole core:
                # One ring is a five-membered aromatic ring with nitrogen
                # The other is a six-membered aromatic ring
                ring_i = rings[i]
                ring_j = rings[j]

                conditions = [
                    (ring_i['size'] == 5 and ring_i['aromatic'] and ring_i['has_nitrogen']
                     and ring_j['size'] == 6 and ring_j['aromatic']),
                    (ring_j['size'] == 5 and ring_j['aromatic'] and ring_j['has_nitrogen']
                     and ring_i['size'] == 6 and ring_i['aromatic'])
                ]
                if any(conditions):
                    found_indole = True
                    break
        if found_indole:
            break

    if found_indole:
        return True, "Contains indole skeleton and nitrogen atoms (indole alkaloid)"
    else:
        return False, "No indole skeleton found"