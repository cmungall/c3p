"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:51294 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is a sphingosine with an acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-acylated nitrogen pattern: [NX3;H0][CX3](=O)[#6]
    n_acyl_pattern = Chem.MolFromSmarts("[NX3;H0][CX3](=O)[#6]")
    n_acyl_matches = mol.GetSubstructMatches(n_acyl_pattern)
    if not n_acyl_matches:
        return False, "No N-acylated nitrogen found"

    # For each acylated nitrogen, check for sphingosine backbone
    for match in n_acyl_matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)

        # Find alpha carbon (C2) connected to nitrogen (excluding carbonyl carbon)
        neighbor_atoms = [a for a in nitrogen_atom.GetNeighbors() if a.GetAtomicNum() == 6 and a.GetIdx() != match[1]]
        if not neighbor_atoms:
            continue
        alpha_carbon = neighbor_atoms[0]

        # Check if alpha carbon has hydroxyl group
        has_alpha_oh = False
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                has_alpha_oh = True
                break
        if not has_alpha_oh:
            continue

        # Find beta carbon (C3) connected to alpha carbon (excluding nitrogen)
        beta_neighbors = [a for a in alpha_carbon.GetNeighbors() if a.GetAtomicNum() == 6 and a.GetIdx() != nitrogen_idx]
        if not beta_neighbors:
            continue
        beta_carbon = beta_neighbors[0]

        # Check if beta carbon has hydroxyl group
        has_beta_oh = False
        for neighbor in beta_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                has_beta_oh = True
                break
        if not has_beta_oh:
            continue

        # From beta carbon, check for long hydrocarbon chain (length ≥12)
        chain_length = 0
        visited = set()
        atoms_to_visit = [beta_carbon]
        while atoms_to_visit:
            current_atom = atoms_to_visit.pop()
            idx = current_atom.GetIdx()
            if idx in visited:
                continue
            visited.add(idx)
            if current_atom.GetAtomicNum() == 6:
                chain_length += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        atoms_to_visit.append(neighbor)
        if chain_length >= 12:
            return True, "Molecule is an N-acylsphingosine with sphingosine backbone"
        else:
            return False, f"Hydrocarbon chain too short (length {chain_length})"

    return False, "Does not match N-acylsphingosine structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51294',
                              'name': 'N-acylsphingosine',
                              'definition': 'The parent compounds of the ceramide family, composed of '
                                            'sphingosine having an unspecified fatty acyl group attached '
                                            'to the nitrogen.',
                              'parents': ['CHEBI:76165', 'CHEBI:76158']},
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
        'num_true_positives': 145,
        'num_false_positives': 5,
        'num_true_negatives': 182406,
        'num_false_negatives': 28,
        'num_negatives': None,
        'precision': 0.9668874172185431,
        'recall': 0.838150289017341,
        'f1': 0.8978675645342312,
        'accuracy': 0.9998195460913245}