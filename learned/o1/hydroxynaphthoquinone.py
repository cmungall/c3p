"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:51793 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is any naphthoquinone where the naphthoquinone moiety is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the naphthoquinone SMARTS pattern (1,4-naphthoquinone core)
    naphthoquinone_smarts = 'O=C1C=CC=CC2=CC=CC(=O)C12'
    naphthoquinone_mol = Chem.MolFromSmarts(naphthoquinone_smarts)

    # Find substructure matches of the naphthoquinone moiety
    matches = mol.GetSubstructMatches(naphthoquinone_mol)
    if not matches:
        return False, "No naphthoquinone moiety found"

    # For each match, check for hydroxy substitutions on the naphthoquinone ring
    for match in matches:
        # Get the atoms involved in the naphthoquinone moiety
        naphthoquinone_atoms = set(match)

        # Identify ketone carbons (carbons double-bonded to oxygen)
        ketone_carbons = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon atom
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_carbons.append(idx)
                        break

        # Identify ring carbons excluding ketone carbons
        ring_carbons = [idx for idx in match if idx not in ketone_carbons and
                        mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and
                        mol.GetAtomWithIdx(idx).IsInRing()]

        # Check for hydroxy groups attached to ring carbons
        has_hydroxy = False
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # Look for oxygen atom with degree 1 (hydroxy group)
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    # Check if oxygen is connected to hydrogen
                    if neighbor.GetTotalNumHs() > 0:
                        has_hydroxy = True
                        break
            if has_hydroxy:
                break

        if has_hydroxy:
            return True, "Contains naphthoquinone moiety substituted with at least one hydroxy group"

    return False, "No hydroxy substitution on naphthoquinone moiety found"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51793',
        'name': 'hydroxynaphthoquinone',
        'definition': 'Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.',
        'parents': ['CHEBI:36124']  # CHEBI:36124 is naphthoquinone
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}