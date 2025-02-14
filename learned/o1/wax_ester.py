"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:50447 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of
    the carboxy group of a fatty acid with the alcoholic hydroxy group of
    a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ester_matches) == 0:
        return False, "No ester functional group found"
    elif len(ester_matches) > 1:
        return False, f"Multiple ester groups found ({len(ester_matches)}), molecule should have only one ester linkage"

    # Get the atoms involved in the ester group
    ester_atoms = ester_matches[0]
    carbonyl_c_idx = ester_atoms[0]
    carbonyl_o_idx = ester_atoms[1]
    alkoxy_o_idx = ester_atoms[2]

    # Initialize sets for atoms in the acyl and alkoxy chains
    acyl_atoms = set()
    alkoxy_atoms = set()

    # Traverse the acyl chain (excluding ester oxygen)
    def traverse_acyl(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
            if n_idx != carbonyl_o_idx and n_idx not in acyl_atoms and bond is not None:
                acyl_atoms.add(n_idx)
                if neighbor.GetAtomicNum() == 6:  # carbon
                    traverse_acyl(n_idx)

    traverse_acyl(carbonyl_c_idx)

    # Traverse the alkoxy chain (excluding carbonyl carbon)
    def traverse_alkoxy(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
            if n_idx != carbonyl_c_idx and n_idx not in alkoxy_atoms and bond is not None:
                alkoxy_atoms.add(n_idx)
                if neighbor.GetAtomicNum() == 6:  # carbon
                    traverse_alkoxy(n_idx)

    traverse_alkoxy(alkoxy_o_idx)

    # Count the number of carbons in each chain
    acyl_carbons = [idx for idx in acyl_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    alkoxy_carbons = [idx for idx in alkoxy_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

    acyl_chain_length = len(acyl_carbons)
    alkoxy_chain_length = len(alkoxy_carbons)

    # Check if both chains are long enough (≥12 carbons)
    if acyl_chain_length < 12:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons), need at least 12"
    if alkoxy_chain_length < 12:
        return False, f"Alkoxy chain too short ({alkoxy_chain_length} carbons), need at least 12"

    # Check for other functional groups
    # Define unwanted functional groups
    unwanted_funcs = {
        'alcohol': '[OX2H]',
        'aldehyde': '[CX3H1](=O)',
        'ketone': '[#6][CX3](=O)[#6]',
        'carboxylic acid': '[CX3](=O)[OX1H1]',
        'amine': '[NX3;H2,H1;!$(NC=O)]',
        'amide': '[NX3][CX3](=O)[#6]',
        'nitrile': '[CX2]#N',
        'nitro': '[NX3](=O)=O',
        'thiol': '[SX2H]',
        'sulfide': '[SX2][#6]',
        'halide': '[F,Cl,Br,I]',
    }

    # Exclude ester group when searching for other functional groups
    ester_atom_indices = set(ester_atoms)

    for name, smarts in unwanted_funcs.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # Check if any of the atoms in the match are not part of the ester group
            if not ester_atom_indices.issuperset(match):
                return False, f"Unwanted functional group ({name}) found"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, which is not characteristic of wax esters"

    return True, "Molecule is a wax ester with long-chain fatty acid and fatty alcohol"