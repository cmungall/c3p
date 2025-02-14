"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:32988 amine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two,
    or three hydrogen atoms with hydrocarbyl groups (alkyl or aryl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to indicate if an amine nitrogen is found
    amine_found = False

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip non-nitrogen atoms

        # Exclude nitrogens that are part of aromatic rings
        if atom.GetIsAromatic():
            continue

        # Exclude amide nitrogens (N connected to C=O)
        is_amide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if the carbon is part of a carbonyl group
                for carb_neighbor in neighbor.GetNeighbors():
                    if carb_neighbor.GetAtomicNum() == 8:  # Oxygen
                        bond = neighbor.GetBondBetweenAtoms(neighbor.GetIdx(), carb_neighbor.GetIdx())
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue  # Skip amide nitrogens

        # Exclude nitro groups (N bonded to two oxygens)
        is_nitro = False
        o_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(o_neighbors) >= 2:
            is_nitro = True
        if is_nitro:
            continue  # Skip nitro nitrogens

        # Check if nitrogen is sp3-hybridized
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue  # Skip if not sp3

        # Check that all substituents (excluding hydrogens) are carbons
        all_carbons = True
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                all_carbons = False
                break
        if not all_carbons:
            continue  # Skip if any substituent is not carbon

        # Nitrogen meets criteria for an amine
        amine_found = True
        break  # No need to check further

    if amine_found:
        return True, "Contains nitrogen atom(s) fitting the definition of an amine"
    else:
        return False, "No nitrogen atoms matching amine criteria found"