"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40), derived from psi,psi-carotene by various modifications
    such as hydrogenation, dehydrogenation, cyclization, oxidation, or combinations thereof.
    Retinoids are excluded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude retinoids by checking for retinoid substructure
    # Retinoids typically have a beta-ionone ring connected to a polyene chain ending with an aldehyde or acid
    retinoid_pattern = Chem.MolFromSmarts("C1=CC=CC=C1CC=CC(=O)[O,N]")
    if mol.HasSubstructMatch(retinoid_pattern):
        return False, "Molecule is a retinoid, which is excluded"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 41:
        return False, f"Molecule has {c_count} carbon atoms, expected between 35 and 41 for carotenoids"

    # Detect long conjugated polyene chain (at least 7 conjugated double bonds)
    # Use a recursive SMARTS pattern to match conjugated double bonds
    polyene_pattern = Chem.MolFromSmarts("[!$([#6]=,#[#6])]~[#6]=[#6]~[#6]=[#6]~[#6]=[#6]~[#6]=[#6]~[#6]=[#6]~[#6]=[#6]~[#6]")
    matches = mol.GetSubstructMatches(polyene_pattern)
    if not matches:
        return False, "No long conjugated polyene chain found (at least 7 conjugated double bonds required)"

    # Consider the molecule a carotenoid
    return True, "Molecule matches criteria for a carotenoid"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids (C40), formally derived from the acyclic parent, psi,psi-carotene by hydrogenation, dehydrogenation, cyclization, oxidation, or combination of these processes. This class includes carotenes, xanthophylls and certain compounds that arise from rearrangement of the skeleton of psi,psi-carotene or by loss of part of this structure. Retinoids are excluded.',
                          'parents': ['CHEBI:26859']},
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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None }