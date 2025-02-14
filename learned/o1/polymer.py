"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:26189 polymer
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    This is done by analyzing larger atom environments to detect repeating units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """
    from collections import defaultdict

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - polymers typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:  # Adjust this threshold as needed
        return False, f"Molecular weight {mol_wt:.2f} Da is too low to be a polymer"

    # Dictionary to count substructures
    substruct_counts = defaultdict(int)
    radius = 4  # Increased radius for larger environments

    for atom_idx in range(mol.GetNumAtoms()):
        # Get the atom environment for the given radius
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
        if not env:
            continue  # Skip if no environment is found
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        # Generate the canonical SMILES of the substructure
        smiles_submol = Chem.MolToSmiles(submol, canonical=True)
        substruct_counts[smiles_submol] += 1

    if not substruct_counts:
        return False, "No substructures found for analysis"

    # Find the maximum count of any substructure
    max_count = max(substruct_counts.values())
    unique_substructs = len(substruct_counts)

    # Thresholds for repeating units and diversity
    repeating_unit_threshold = 3  # Number of times a substructure should repeat
    diversity_threshold = 0.5  # Ratio of total substructures to unique substructures

    # Calculate diversity ratio
    total_substructs = sum(substruct_counts.values())
    diversity_ratio = unique_substructs / total_substructs

    if max_count >= repeating_unit_threshold and diversity_ratio < diversity_threshold:
        return True, f"Molecule has repeating substructures occurring {max_count} times with low diversity ({diversity_ratio:.2f}), likely a polymer"
    else:
        return False, "No significant repeating units found or high structural diversity"