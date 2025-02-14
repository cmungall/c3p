"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29232 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over sulfur atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check if sulfur is bonded to at least one carbon atom
            bonded_to_carbon = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    bonded_to_carbon = True
                    break
            if not bonded_to_carbon:
                continue  # Not a thiol if sulfur is not bonded to carbon

            # Get total number of hydrogens (implicit and explicit) on sulfur
            total_hydrogens = atom.GetTotalNumHs()
            if total_hydrogens != 1:
                continue  # Thiol sulfur must have exactly one hydrogen

            # Check that sulfur is not double-bonded (no multiple bonds)
            multiple_bonds = False
            for bond in atom.GetBonds():
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    multiple_bonds = True
                    break
            if multiple_bonds:
                continue  # Not a thiol if sulfur has multiple bonds

            # Check that sulfur is not part of a disulfide bond (S-S)
            bonded_to_sulfur = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 16:  # Sulfur atom
                    bonded_to_sulfur = True
                    break
            if bonded_to_sulfur:
                continue  # Not a thiol if sulfur is bonded to another sulfur

            # All conditions met; molecule contains a thiol group
            return True, "Contains thiol group (-SH) attached to a carbon atom"

    return False, "No thiol group (-SH) attached to a carbon atom found"