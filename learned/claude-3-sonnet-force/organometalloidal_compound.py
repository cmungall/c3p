"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: CHEBI:33709 organometalloidal compound

An organometalloidal compound is defined as a compound having bonds between one or more metalloid
atoms and one or more carbon atoms of an organyl group.
"""

from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get list of metalloid atomic numbers
    metalloids = [5, 14, 32, 33, 51, 52, 83, 84]

    # Check if molecule contains any metalloid atoms
    has_metalloid = any(atom.GetAtomicNum() in metalloids for atom in mol.GetAtoms())
    if not has_metalloid:
        return False, "No metalloid atoms found"

    # Check if any metalloid atoms are bonded to carbon atoms
    for bond in mol.GetBonds():
        atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if (atom1.GetAtomicNum() in metalloids and atom2.GetAtomicNum() == 6) or \
           (atom2.GetAtomicNum() in metalloids and atom1.GetAtomicNum() == 6):
            return True, "Contains bonds between metalloid and carbon atoms"

    return False, "No bonds between metalloid and carbon atoms found"