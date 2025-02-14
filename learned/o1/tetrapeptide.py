"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule containing four amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the peptide bond pattern (backbone peptide bonds)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    peptide_bond_atoms = []
    for match in peptide_bond_matches:
        c_idx, n_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Check if N atom is connected to a C-alpha carbon (sp3 carbon with one hydrogen)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetHybridization() == Chem.HybridizationType.SP3:
                num_hydrogens = neighbor.GetTotalNumHs()
                if num_hydrogens == 1:
                    # This is likely a C-alpha carbon
                    peptide_bond_atoms.append((c_idx, n_idx))
                    break

    num_peptide_bonds = len(peptide_bond_atoms)

    if num_peptide_bonds != 3:
        return False, f"Found {num_peptide_bonds} peptide bonds, expected 3 for tetrapeptide"

    # The number of residues is peptide bonds + 1
    num_residues = num_peptide_bonds + 1

    if num_residues != 4:
        return False, f"Found {num_residues} amino acid residues, expected 4 for tetrapeptide"

    return True, "Molecule is a tetrapeptide with 4 amino acid residues connected via peptide bonds"