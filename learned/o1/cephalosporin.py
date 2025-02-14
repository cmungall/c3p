"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:23066 cephalosporin
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are beta-lactam antibiotics characterized by a beta-lactam ring fused to a six-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Lists to store ring indices
    beta_lactam_rings = []
    dihydrothiazine_rings = []

    # Identify beta-lactam and dihydrothiazine rings
    for ring in atom_rings:
        if len(ring) == 4:
            # Potential beta-lactam ring
            num_N = 0
            has_carbonyl = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    num_N += 1
            if num_N == 1:
                # Check for carbonyl group (C=O) in ring
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    has_carbonyl = True
                if has_carbonyl:
                    beta_lactam_rings.append(set(ring))
        elif len(ring) == 6:
            # Potential dihydrothiazine ring
            num_N = 0
            num_S = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    num_N += 1
                elif atom.GetAtomicNum() == 16:
                    num_S += 1
            if num_N == 1 and num_S == 1:
                dihydrothiazine_rings.append(set(ring))

    # Check for fused rings
    for beta_ring in beta_lactam_rings:
        for dihydro_ring in dihydrothiazine_rings:
            shared_atoms = beta_ring & dihydro_ring
            if len(shared_atoms) >= 2:
                # Rings are fused
                return True, "Contains fused beta-lactam and dihydrothiazine rings characteristic of cephalosporins"

    return False, "Cephalosporin core structure not found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring. Although cephalosporins are among the most commonly used antibiotics in the treatment of routine infections, and their use is increasing over time, they can cause a range of hypersensitivity reactions, from mild, delayed-onset cutaneous reactions to life-threatening anaphylaxis in patients with immunoglobulin E (IgE)-mediated allergy.',
    },
}