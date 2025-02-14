"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amine nitrogens (excluding amides, nitriles, etc.)
    amine_nitrogens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Exclude nitrogen atoms that are double or triple bonded
            if atom.GetDegree() > 0 and atom.GetHybridization() == Chem.HybridizationType.SP3:
                # Exclude amide nitrogens (connected to carbonyl carbon)
                is_amide = False
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 6:
                        for nbr_bond in neighbor.GetBonds():
                            if nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr_bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                                is_amide = True
                if not is_amide:
                    amine_nitrogens.append(atom)

    if not amine_nitrogens:
        return False, "No amine group found"

    # For each amine nitrogen, search for an aromatic ring connected via an alkyl chain
    for nitrogen in amine_nitrogens:
        visited = set()
        queue = [(nitrogen.GetIdx(), 0)]  # tuple of (atom idx, distance)
        max_distance = 6  # Limit the search to avoid infinite loops
        found_aromatic = False
        while queue:
            current_idx, distance = queue.pop(0)
            if current_idx in visited or distance > max_distance:
                continue
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            if current_atom.GetIdx() != nitrogen.GetIdx():
                # Skip if the atom is another heteroatom
                if current_atom.GetAtomicNum() != 6:
                    continue
                # Skip if the atom is not sp3 hybridized (to stay within alkyl chain)
                if current_atom.GetHybridization() != Chem.HybridizationType.SP3:
                    continue
            for bond in current_atom.GetBonds():
                neighbor = bond.GetOtherAtom(current_atom)
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in visited:
                    continue
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                if neighbor.GetAtomicNum() == 6:
                    if neighbor.GetIsAromatic():
                        found_aromatic = True
                        break
                    else:
                        queue.append((neighbor_idx, distance + 1))
                elif neighbor.GetAtomicNum() in [7,8,16]:
                    continue  # Skip other heteroatoms
            if found_aromatic:
                break
        if found_aromatic:
            return True, "Amine group connected via alkyl chain to an aromatic ring"

    return False, "No aralkylamine substructure found"