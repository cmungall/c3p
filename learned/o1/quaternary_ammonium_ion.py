"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35279 quaternary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a nitrogen atom with a +1 charge, bonded to four univalent
    (usually organyl) groups (carbon atoms), and not part of any ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetFormalCharge() == 1:  # Formal charge +1
                if atom.GetDegree() == 4:  # Bonded to four atoms
                    if not atom.IsInRing():  # Not in a ring
                        neighbors = atom.GetNeighbors()
                        all_single_bonds = True
                        all_carbon_neighbors = True
                        for bond in atom.GetBonds():
                            if bond.GetBondType() != rdchem.BondType.SINGLE:
                                all_single_bonds = False
                                break
                        for neighbor in neighbors:
                            if neighbor.GetAtomicNum() != 6:  # Carbon atom
                                all_carbon_neighbors = False
                                break
                        if all_single_bonds and all_carbon_neighbors:
                            return True, "Contains a quaternary ammonium ion (nitrogen with +1 charge bonded to four carbon atoms via single bonds)"
    return False, "No quaternary ammonium ion found in the molecule"