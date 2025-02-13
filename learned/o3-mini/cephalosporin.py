"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: Cephalosporin
A cephalosporin is characterized by a cephem nucleus: a beta-lactam (4-membered cyclic amide) fused to a 6-membered dihydrothiazine ring (which contains a sulfur). This program attempts to detect that fused bicyclic system.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    The algorithm looks for a 4-membered beta-lactam ring (with an amide nitrogen and a carbonyl carbon)
    fused to a 6-membered ring containing at least one sulfur atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a cephalosporin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or not ring_info.AtomRings():
        return False, "No ring information could be detected"

    atom_rings = ring_info.AtomRings()
    
    # Find candidate beta-lactam rings: 4-membered rings with one nitrogen and one carbonyl carbon (C double-bonded to O)
    beta_lactam_rings = []
    for ring in atom_rings:
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Count nitrogen atoms in the ring
            n_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            # Check for a carbonyl group: find a carbon atom in the ring which is bonded via a double bond to an oxygen
            has_carbonyl = False
            for atom in atoms_in_ring:
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Check if bond is double and neighbor is oxygen
                        if bond.GetBondTypeAsDouble() == 2.0:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            # We expect at least one nitrogen and one carbonyl carbon in the beta-lactam ring
            if n_nitrogen >= 1 and has_carbonyl:
                beta_lactam_rings.append(set(ring))
    
    if not beta_lactam_rings:
        return False, "No beta-lactam (4-membered amide) ring detected"
    
    # Find candidate 6-membered rings that contain at least one sulfur atom (typical for dihydrothiazine)
    six_membered_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if any(atom.GetAtomicNum() == 16 for atom in atoms_in_ring):
                six_membered_rings.append(set(ring))
    
    if not six_membered_rings:
        return False, "No 6-membered ring containing sulfur detected (dihydrothiazine ring missing)"
    
    # Check for fused rings: the beta-lactam ring and the 6-membered ring must share at least 2 atoms
    for beta_ring in beta_lactam_rings:
        for six_ring in six_membered_rings:
            common_atoms = beta_ring.intersection(six_ring)
            if len(common_atoms) >= 2:
                # Optionally, one may add further checks such as verifying a carboxylate group is present.
                return True, "Molecule has a beta-lactam ring fused with a 6-membered sulfur-containing ring (cephem nucleus)"
    
    return False, "No fused beta-lactam and 6-membered sulfur-containing ring system detected"

# Example usage:
if __name__ == '__main__':
    # An example SMILES from the provided list of cephalosporins
    test_smiles = "CC1=C(C(=NN1CC(C)C(=O)NC2C3N(C2=O)C(=C(CS3)COC(=O)C)C(=O)O)C)Cl"
    result, reason = is_cephalosporin(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)