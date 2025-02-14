"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a compound where a nitrogen atom is sp3-hybridized,
    bonded to three hydrocarbyl groups (carbon atoms), not part of an amide, imine,
    nitrile, or aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to indicate if a tertiary amine nitrogen is found
    is_tertiary_amine = False

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Exclude charged nitrogen atoms
            if atom.GetFormalCharge() != 0:
                continue

            # Exclude aromatic nitrogen atoms
            if atom.GetIsAromatic():
                continue

            # Check hybridization (should be sp3)
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                continue

            # Nitrogen should have degree 3 (three neighbors)
            if atom.GetDegree() != 3:
                continue

            # Nitrogen should have no hydrogens
            if atom.GetTotalNumHs(includeNeighbors=True) != 0:
                continue

            # Check that all neighbors are carbon atoms and not involved in multiple bonds
            is_valid = True
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() != 6:  # Not carbon
                    is_valid = False
                    break
                # Exclude multiple bonds
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    is_valid = False
                    break
                # Exclude if neighbor is part of a carbonyl group (C=O) or imine (C=N)
                for nbr_bond in neighbor.GetBonds():
                    if nbr_bond.GetBondType() == rdchem.BondType.DOUBLE:
                        bonded_atom = nbr_bond.GetOtherAtom(neighbor)
                        if bonded_atom.GetAtomicNum() in [7, 8]:  # Nitrogen or Oxygen
                            is_valid = False
                            break
                if not is_valid:
                    break

            if is_valid:
                is_tertiary_amine = True
                break  # No need to check further

    if is_tertiary_amine:
        return True, "Contains a sp3-hybridized nitrogen atom bonded to three carbon atoms (tertiary amine)"
    else:
        return False, "No tertiary amine group found"