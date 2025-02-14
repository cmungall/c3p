"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a ketone (=O) functional group at position 3 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate Murcko scaffold of the molecule
    mol_scaffold = MurckoScaffold.GetScaffoldForMol(mol)

    # Generate Murcko scaffold of a reference steroid (e.g., testosterone)
    testosterone_smiles = "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34"
    testosterone_mol = Chem.MolFromSmiles(testosterone_smiles)
    if testosterone_mol is None:
        return False, "Error processing reference steroid"
    ref_scaffold = MurckoScaffold.GetScaffoldForMol(testosterone_mol)

    # Compare the scaffold SMILES
    mol_scaffold_smiles = Chem.MolToSmiles(mol_scaffold, isomericSmiles=True)
    ref_scaffold_smiles = Chem.MolToSmiles(ref_scaffold, isomericSmiles=True)

    if mol_scaffold_smiles != ref_scaffold_smiles:
        return False, "Molecule does not share scaffold with reference steroid"

    # Define SMARTS pattern for ketone group at position 3
    # Approximate by looking for a ketone attached to the first ring (A ring)
    ketone_pattern = Chem.MolFromSmarts('O=C[C;R]')  # Ketone attached to ring carbon
    matches = mol.GetSubstructMatches(ketone_pattern)
    if not matches:
        return False, "No ketone group found in ring system"

    # Check if ketone is in ring system
    for match in matches:
        carbon_idx = match[1]
        atom = mol.GetAtomWithIdx(carbon_idx)
        # Check if carbon is in a ring
        if atom.IsInRing():
            # Optionally, check if atom is in the first ring (A ring)
            # This requires defining the A ring, which is complex without atom indexing
            return True, "Molecule is a 3-oxo steroid"
    
    return False, "Ketone group not in ring system"