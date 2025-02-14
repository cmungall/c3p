"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: Vitamin D compounds
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are secosteroids with a broken B-ring, a conjugated triene system,
    hydroxyl groups at specific positions, and a hydrocarbon side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for secosteroid core (rings C and D fused together)
    # Vitamin D compounds have fused cyclohexane (ring C) and cyclopentane (ring D)
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()
    ring_sizes = [len(ring) for ring in atom_rings]

    # Look for one ring of size 6 (ring C) and one ring of size 5 (ring D)
    if ring_sizes.count(5) < 1 or ring_sizes.count(6) < 1:
        return False, "Secosteroid core not found (missing fused rings C and D)"

    # Check if rings C and D are fused (share two atoms)
    fused = False
    for i, ring1 in enumerate(atom_rings):
        for j, ring2 in enumerate(atom_rings):
            if i >= j:
                continue
            if len(ring1) == 6 and len(ring2) == 5:
                shared_atoms = set(ring1) & set(ring2)
                if len(shared_atoms) >= 2:
                    fused = True
                    break
    if not fused:
        return False, "Rings C and D are not fused"

    # Check for conjugated triene system in the open A-ring
    triene_smarts = 'C=C-C=C-C=C'  # Conjugated triene system
    triene_pattern = Chem.MolFromSmarts(triene_smarts)
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system found"

    # Check for at least two hydroxyl groups
    hydroxyl_smarts = '[OX2H]'  # Hydroxyl group
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts(hydroxyl_smarts)))
    if hydroxyl_count < 2:
        return False, f"Found {hydroxyl_count} hydroxyl group(s), need at least 2"

    # Estimate molecular weight for a vitamin D range (approximately 380-450 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} Da not in typical range for vitamin D"

    return True, "Molecule matches key structural features of vitamin D"