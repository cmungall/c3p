"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoethanolamine is a glycerophosphoethanolamine with an acyl group 
    at the sn-1 position and (R)-configuration at the chiral center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to correctly perceive chiral centers
    mol = Chem.AddHs(mol)
    
    # Find chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found"

    # Identify sn-2 carbon (should be a carbon attached to two oxygens and one carbon)
    sn2_atom_idx = None
    for idx, config in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue  # Not a carbon
        neighbors = atom.GetNeighbors()
        o_count = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 8)
        c_count = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
        if o_count == 2 and c_count == 1:
            sn2_atom_idx = idx
            sn2_config = config
            break
    if sn2_atom_idx is None:
        return False, "No appropriate chiral center found for sn-2 position"

    # Check sn-2 configuration
    if sn2_config != 'R':
        return False, f"Chiral center at sn-2 position is not (R)-configuration (found {sn2_config})"

    # Check for ester linkage at sn-1 position
    # sn-1 position is primary alcohol connected to acyl group via ester linkage
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found at sn-1 position"

    # Verify that the O in ester linkage is connected to sn-1 carbon (neighbor of sn-2 carbon)
    sn1_found = False
    for match in ester_matches:
        ester_o_idx = match[0]  # Index of oxygen in ester linkage
        ester_o = mol.GetAtomWithIdx(ester_o_idx)
        for neighbor in ester_o.GetNeighbors():
            if neighbor.GetIdx() == sn2_atom_idx:
                sn1_found = True
                break
        if sn1_found:
            break
    if not sn1_found:
        return False, "Ester linkage at sn-1 position not connected properly"

    # Check for phosphoethanolamine group at sn-3 position
    # sn-3 position is primary alcohol connected to phosphoethanolamine group
    phospho_pattern = Chem.MolFromSmarts("O[P](=O)(OCCN)O")
    phospho_matches = mol.GetSubstructMatches(phospho_pattern)
    if not phospho_matches:
        return False, "No phosphoethanolamine group found at sn-3 position"

    # Verify that the O in phosphoester linkage is connected to sn-3 carbon (neighbor of sn-2 carbon)
    sn3_found = False
    for match in phospho_matches:
        phospho_o_idx = match[3]  # Index of oxygen connecting phospho group to glycerol backbone
        phospho_o = mol.GetAtomWithIdx(phospho_o_idx)
        for neighbor in phospho_o.GetNeighbors():
            if neighbor.GetIdx() == sn2_atom_idx:
                sn3_found = True
                break
        if sn3_found:
            break
    if not sn3_found:
        return False, "Phosphoethanolamine group at sn-3 position not connected properly"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine with (R)-configuration"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '1-acyl-sn-glycero-3-phosphoethanolamine',
        'definition': 'A 1-O-acylglycerophosphoethanolamine having (R)-configuration.',
        'parents': []
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
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}