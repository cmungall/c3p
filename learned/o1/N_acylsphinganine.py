"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:52634 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is a ceramide consisting of sphinganine where one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sphinganine backbone pattern:
    # Long-chain amino diol with hydroxyls at positions 1 and 3, amino group at position 2
    sphinganine_smarts = """
    [C;!$(C=*)]            # Terminal carbon (position 1)
    [C@H](O)               # Chiral carbon with hydroxyl (position 1)
    [C@H](N)               # Chiral carbon with amino group (position 2)
    [C@H](O)               # Chiral carbon with hydroxyl (position 3)
    [CH2]                  # Remaining methylene carbons in the chain
    {[CH2],[CH]}*          # Allowing for variable chain lengths
    [CH3]                  # Terminal methyl group
    """
    sphinganine_pattern = Chem.MolFromSmarts(sphinganine_smarts)
    if sphinganine_pattern is None:
        return False, "Invalid sphinganine SMARTS pattern"

    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"

    # Define N-acyl amide group pattern attached to the amino group
    n_acyl_smarts = "[N][C](=O)[C;R0]([C;R0])([C;R0])[C;R0]"  # Amide with aliphatic chain
    n_acyl_pattern = Chem.MolFromSmarts(n_acyl_smarts)
    if n_acyl_pattern is None:
        return False, "Invalid N-acyl SMARTS pattern"

    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group (amide linkage) found attached to the amino nitrogen"

    # Check that the acyl chain is long (fatty acyl group)
    acyl_chain_smarts = "C(=O)[C;R0].[C;R0]"  # Carbonyl carbon attached to a long chain
    acyl_chain_pattern = Chem.MolFromSmarts(acyl_chain_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)

    # Estimate acyl chain length
    acyl_lengths = []
    for match in acyl_matches:
        carbonyl_c_idx = match[0]
        acyl_c_idx = match[1]
        # Traverse the acyl chain to count carbons
        visited = set()
        queue = [acyl_c_idx]
        length = 0
        while queue:
            atom_idx = queue.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                length += 1
                for neighbor in atom.GetNeighbors():
                    nbr_idx = neighbor.GetIdx()
                    if neighbor.GetAtomicNum() == 6 and nbr_idx != carbonyl_c_idx:
                        queue.append(nbr_idx)
        acyl_lengths.append(length)

    if not acyl_lengths or max(acyl_lengths) < 12:
        return False, "Acyl chain is too short to be a fatty acyl group"

    # Confirm that the N-acyl group is connected to the amino nitrogen of sphinganine
    # Find the amino nitrogen atom in sphinganine backbone
    amino_nitrogen_smarts = "[C@H](N)[C@H](O)"
    amino_nitrogen_pattern = Chem.MolFromSmarts(amino_nitrogen_smarts)
    amino_matches = mol.GetSubstructMatches(amino_nitrogen_pattern)
    if not amino_matches:
        return False, "No amino nitrogen found in sphinganine backbone"

    # Ensure the amino nitrogen is part of the amide linkage
    amide_nitrogen_smarts = "[N][C](=O)"
    amide_nitrogen_pattern = Chem.MolFromSmarts(amide_nitrogen_smarts)
    if not mol.HasSubstructMatch(amide_nitrogen_pattern):
        return False, "Amino nitrogen is not part of an amide linkage"

    # Check total chain length (sphinganine + fatty acyl chain)
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() <= 4:
            total_carbons += 1
    if total_carbons < 18:
        return False, "Total carbon count is too low for N-acylsphinganine"

    return True, "Contains sphinganine backbone with N-acyl group forming an amide bond"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52634',
        'name': 'N-acylsphinganine',
        'definition': "A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.",
        'parents': ['CHEBI:33598', 'CHEBI:33709']
    },
    'config': {}
}