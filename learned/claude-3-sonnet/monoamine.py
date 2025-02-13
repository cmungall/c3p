"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI:35555 monoamine

A monoamine is defined as an aralylamino compound which contains one amino group connected to an aromatic ring
by a two-carbon chain. Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan,
and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for an aromatic ring
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if mol.GetRingAtoms(ring) == list(aromatic_atoms)]
    if not aromatic_rings:
        return False, "No aromatic ring found"
    
    # Look for an amino group (-NH2 or -NH3+)
    amino_pattern = Chem.MolFromSmarts("[N;H2,H3]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"
    
    # Check if the amino group is connected to the aromatic ring by a two-carbon chain
    for amino_idx in amino_matches:
        chain_atoms = mol.GetAtomWithIdx(amino_idx).FindAllPathsOfLengthN(aromatic_atoms, 3, usehsenv=False)
        if chain_atoms:
            chain_atoms = [mol.GetAtomWithIdx(idx) for paths in chain_atoms for idx in paths[1:-1]]
            if all(atom.GetAtomicNum() == 6 for atom in chain_atoms):
                return True, "Contains an amino group connected to an aromatic ring by a two-carbon chain"
    
    return False, "Amino group not connected to aromatic ring by a two-carbon chain"