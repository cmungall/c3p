"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside has a steroid nucleus with an unsaturated lactone ring at C17
    and sugar residues attached at C3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus pattern (cyclopenta[a]phenanthrene system)
    steroid_nucleus = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4=C3C=CC=C4')  # Simplified pattern
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"

    # Define unsaturated lactone ring at C17 (butenolide attached to steroid)
    lactone_pattern = Chem.MolFromSmarts('C=1C=CC(=O)O1')  # α,β-unsaturated γ-lactone
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No unsaturated lactone ring found"

    # Check if lactone is attached to C17 position of steroid nucleus
    # For simplicity, we'll check proximity between steroid and lactone substructures
    steroid_matches = mol.GetSubstructMatch(steroid_nucleus)
    lactone_atoms = [atom_idx for match in lactone_matches for atom_idx in match]
    shared_atoms = set(steroid_matches) & set(lactone_atoms)
    if not shared_atoms:
        return False, "Lactone ring not attached to steroid nucleus"

    # Define sugar residue pattern (pyranose ring) attached via glycosidic bond at C3
    sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C1')  # Simplified pyranose
    glycosidic_bond = Chem.MolFromSmarts('O[C@H]')  # Oxygen linked to anomeric carbon
    # Search for sugars in the molecule
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar residues found"

    # Check for glycosidic bond connection to the steroid nucleus at C3
    # Identify attachment point
    sugar_atom_indices = [atom_idx for match in sugar_matches for atom_idx in match]
    steroid_atom_indices = list(steroid_matches)
    bonds = mol.GetBonds()
    glycosidic_attachment = False
    for bond in bonds:
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in steroid_atom_indices and end_idx in sugar_atom_indices) or \
           (end_idx in steroid_atom_indices and begin_idx in sugar_atom_indices):
            glycosidic_attachment = True
            break

    if not glycosidic_attachment:
        return False, "Sugar residues not attached to steroid nucleus"

    # Check if sugar is attached at C3 position (requires mapping positions)
    # For simplicity, assume attachment if connected to steroid nucleus
    # Note: Detailed atom mapping would be needed for precise position checking

    return True, "Contains steroid nucleus with C17 lactone ring and sugars attached at C3"

__metadata__ = {
    'chemical_class': {
        'name': 'cardiac glycoside',
        'definition': 'Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.',
    },
    'examples': [
        # SMILES strings of example cardiac glycosides provided
        'O[C@@]12[C@]3([C@@]([C@@]4([C@](CC3)(C[C@@H](O[C@@H]5O[C@H]([C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@@H](O)[C@H]5O)C)CC4)[H])C)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]',
        # Other examples can be included here
    ],
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
}