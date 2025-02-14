"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33561 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation
    of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Sugar ring: five-membered ring with one oxygen and four carbons
    sugar_pattern = Chem.MolFromSmarts("[C;R][O;R][C;R][C;R][C;R]")  # Five-membered ring with one O
    sugar_ring_query = Chem.MolFromSmarts("[$([C;R]1-[O;R]-[C;R]-[C;R]-[C;R]-1)]")  # Ring closure

    # Nitrogenous bases (purines and pyrimidines, including modified bases)
    purine_pattern = Chem.MolFromSmarts("c1ncnc2ncnn12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ccnc(=O)[nH]1")
    base_patterns = [purine_pattern, pyrimidine_pattern]

    # Phosphate group (allowing for various protonation states)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")

    # Identify sugar rings
    sugar_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if len(ring_atoms) == 5:
            o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if o_count == 1 and c_count == 4:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No sugar ring found"

    # Identify nitrogenous bases
    base_found = False
    for base_pattern in base_patterns:
        base_matches = mol.GetSubstructMatches(base_pattern)
        if base_matches:
            base_found = True
            base_atoms = set()
            for match in base_matches:
                base_atoms.update(match)
            break
    if not base_found:
        return False, "No nitrogenous base found"

    # Identify phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    phosphate_atoms = set()
    for match in phosphate_matches:
        phosphate_atoms.update(match)

    # Check connections between sugar and base
    glycosidic_bond_found = False
    for sugar_ring in sugar_rings:
        for carbon_idx in sugar_ring:
            atom = mol.GetAtomWithIdx(carbon_idx)
            if atom.GetAtomicNum() == 6:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7 and neighbor.GetIdx() in base_atoms:
                        # Found glycosidic bond
                        glycosidic_bond_found = True
                        break
                if glycosidic_bond_found:
                    break
        if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
        return False, "No glycosidic bond between sugar and base found"

    # Check connections between sugar and phosphate group
    phosphate_connected_to_sugar = False
    for sugar_ring in sugar_rings:
        for sugar_atom_idx in sugar_ring:
            sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
            if sugar_atom.GetAtomicNum() == 6 or sugar_atom.GetAtomicNum() == 8:
                for neighbor in sugar_atom.GetNeighbors():
                    if neighbor.GetIdx() in phosphate_atoms:
                        # Found connection to phosphate
                        phosphate_connected_to_sugar = True
                        break
                if phosphate_connected_to_sugar:
                    break
        if phosphate_connected_to_sugar:
            break
    if not phosphate_connected_to_sugar:
        return False, "Phosphate group not connected to sugar"

    return True, "Molecule is a nucleotide"