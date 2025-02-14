"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is an N-glycosyl compound that has both a nucleobase
    (normally adenine, guanine, xanthine, thymine, cytosine, or uracil)
    and either a ribose or deoxyribose sugar as functional parents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()

    # Identify sugar rings (five-membered rings containing oxygen)
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            o_count = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == 'O')
            if o_count == 1:
                sugar_rings.append(set(ring))

    if not sugar_rings:
        return False, "No sugar ring found"

    # Identify nucleobase rings (rings containing nitrogen atoms)
    nucleobase_rings = []
    for ring in ring_info.AtomRings():
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        n_count = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == 'N')
        c_count = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == 'C')
        if n_count >= 2 and c_count >= 3:
            nucleobase_rings.append(set(ring))

    if not nucleobase_rings:
        return False, "No nucleobase ring found"

    # Check for N-glycosidic bond between sugar ring and nucleobase ring
    has_nglycosidic_bond = False
    for sugar_ring in sugar_rings:
        for nucleobase_ring in nucleobase_rings:
            # Find bonds between sugar ring and nucleobase ring
            for sugar_atom_idx in sugar_ring:
                sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
                for bond in sugar_atom.GetBonds():
                    neighbor = bond.GetOtherAtom(sugar_atom)
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx in nucleobase_ring and neighbor.GetSymbol() == 'N':
                        # Check for bond between sugar carbon and nucleobase nitrogen
                        if sugar_atom.GetSymbol() == 'C' and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            has_nglycosidic_bond = True
                            break
                if has_nglycosidic_bond:
                    break
            if has_nglycosidic_bond:
                break
        if has_nglycosidic_bond:
            break

    if not has_nglycosidic_bond:
        return False, "No N-glycosidic bond between sugar and nucleobase found"

    return True, "Contains nucleobase attached to sugar via N-glycosidic bond"