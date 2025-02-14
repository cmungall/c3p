"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound is defined as 'A compound having bonds between one or more metalloid atoms and one or more carbon atoms of an organyl group.'

    Metalloid elements considered here are arsenic (As) and antimony (Sb), as they are commonly involved in organometalloidal compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # List of metalloid atomic numbers to consider
    metalloid_atomic_nums = [33, 51]  # As, Sb

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to check for metalloid-carbon bond
    has_metalloid_carbon_bond = False

    # Iterate over bonds
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atomic_num1 = atom1.GetAtomicNum()
        atomic_num2 = atom2.GetAtomicNum()

        # Check if bond is between metalloid and carbon
        if ((atomic_num1 in metalloid_atomic_nums and atomic_num2 == 6) or
            (atomic_num2 in metalloid_atomic_nums and atomic_num1 == 6)):

            # Check if carbon is part of an organyl group (exclude carbonyl carbons)
            if atom1.GetAtomicNum() == 6:
                carbon_atom = atom1
            else:
                carbon_atom = atom2

            # Exclude carbons double-bonded to oxygen (carbonyl groups)
            is_carbonyl = False
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and carbon_atom.GetBondBetweenAtoms(carbon_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    is_carbonyl = True
                    break

            if not is_carbonyl:
                return True, f"Metalloid atom ({atom1.GetSymbol() if atomic_num1 in metalloid_atomic_nums else atom2.GetSymbol()}) bonded to carbon atom of organyl group"

    # If no metalloid-carbon bonds found
    return False, "No metalloid-carbon bonds to organyl groups found"