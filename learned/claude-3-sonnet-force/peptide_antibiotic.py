"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: CHEBI:36352 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are chemically diverse peptides that exhibit antimicrobial properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for peptide backbone
    peptide_pattern = Chem.MolFromSmarts("[NX3]([CH2])[CH2][CX3](=[OX1])[NX3][CH2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide backbone found"
    
    # Check for antimicrobial motifs
    motifs = (
        Chem.MolFromSmarts("[NX3+](C)(C)[CH2][CH2][CH2][NH3+]"),  # polymyxin-like motif
        Chem.MolFromSmarts("[NX3][CH2][CH2]=[CH]C(=[OX1])[NX3]"),  # dehydro amino acid
        Chem.MolFromSmarts("[OX2][CH2][CH2]Cl"),  # chlorinated side chain
        Chem.MolFromSmarts("[OX2][CH2][CH2][NX3+](C)(C)"),  # quaternary amine
        Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][nX3]1"),  # heterocyclic ring
    )
    has_motif = any(mol.HasSubstructMatch(motif) for motif in motifs)
    if not has_motif:
        return False, "No antimicrobial motifs found"
    
    # Check molecular weight - peptide antibiotics typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for peptide antibiotic"
    
    # Count nitrogen, oxygen, and sulfur atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    if n_count < 5 or o_count < 5 or s_count == 0:
        return False, "Atom counts not typical for peptide antibiotic"
    
    return True, "Contains peptide backbone and antimicrobial motifs"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36352',
                          'name': 'peptide antibiotic',
                          'definition': 'A chemically diverse class of peptides '
                                        'that exhibit antimicrobial properties.',
                          'parents': ['CHEBI:35434']},
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
    'num_true_positives': 233,
    'num_false_positives': 24,
    'num_true_negatives': 185924,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9067796610169491,
    'recall': 1.0,
    'f1': 0.9511719511719512,
    'accuracy': 0.999989872103401}