"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for sugar rings (pyranose and furanose)
    pyranose_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1")
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1")

    # Check if molecule contains a sugar ring
    has_sugar_ring = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)
    if not has_sugar_ring:
        return False, "No sugar ring (pyranose or furanose) found"

    # Define pattern for oxygen atoms in ring
    o_in_ring_query = rdqueries.AtomNumEqualsQueryAtom(8)
    o_in_ring_query.SetQueryAtomProp("inRing", 1)

    # Define pattern for sulfur atoms in ring
    s_in_ring_query = rdqueries.AtomNumEqualsQueryAtom(16)
    s_in_ring_query.SetQueryAtomProp("inRing", 1)

    # Count oxygen atoms in ring
    o_in_ring = 0
    s_in_ring = 0
    for atom in mol.GetAtoms():
        if atom.IsInRing():
            if atom.GetAtomicNum() == 8:
                o_in_ring += 1
            elif atom.GetAtomicNum() == 16:
                s_in_ring += 1

    # Check if any ring oxygen is replaced by sulfur
    if s_in_ring == 0:
        ring_heteroatom_replaced = False
    else:
        ring_heteroatom_replaced = True

    # Define pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    # Define pattern for thiol groups (-SH)
    thiol_pattern = Chem.MolFromSmarts("[SX2H]")
    # Define pattern for thioether groups (-SR)
    thioether_pattern = Chem.MolFromSmarts("[SX2]([#6])[#6]")

    # Find hydroxyl groups attached to sugar ring
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    # Find thiol groups attached to sugar ring
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    # Find thioether groups attached to sugar ring
    thioether_matches = mol.GetSubstructMatches(thioether_pattern)

    # Check if any hydroxyl group is replaced by sulfur-containing group
    sulfur_substituted_hydroxyl = len(thiol_matches) + len(thioether_matches) > 0

    # Determine if molecule is a thiosugar
    if ring_heteroatom_replaced or sulfur_substituted_hydroxyl:
        return True, "Molecule is a thiosugar with sulfur substitutions in place of oxygen"
    else:
        return False, "No sulfur substitutions found in place of oxygen in sugar moiety"

__metadata__ = {   'chemical_class': {   'name': 'thiosugar',
                                         'definition': 'A carbohydrate derivative in which one or more of the oxygens or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.'},
                   'config': {   'llm_model_name': 'lbl/claude-sonnet',
                                 'f1_threshold': 0.8,
                                 'max_attempts': 5,
                                 'max_positive_instances': None,
                                 'max_positive_to_test': None,
                                 'max_negative_to_test': None,
                                 'max_positive_in_prompt': 50,
                                 'max_negative_in_prompt': 20,
                                 'max_instances_in_prompt': 100,
                                 'test_proportion': 0.1},
                   'message': None,
                   'attempt': 0,
                   'success': True,
                   'best': True,
                   'error': '',
                   'stdout': None,
                   'num_true_positives': None,
                   'num_false_positives': None,
                   'num_true_negatives': None,
                   'num_false_negatives': None,
                   'num_negatives': None,
                   'precision': None,
                   'recall': None,
                   'f1': None,
                   'accuracy': None}