"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is defined as a ceramide consisting of sphinganine
    in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for sphinganine backbone: NC(O)C(OH)
    sphinganine_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](CO)N")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Check for N-acyl substitution (C(=O) bound to nitrogen)
    n_acyl_pattern = Chem.MolFromSmarts("N[C@H]C=O")
    acyl_matches = mol.GetSubstructMatches(n_acyl_pattern)
    if len(acyl_matches) == 0:
        return False, "No N-acyl substitution found"

    # Check for long carbon chain typically found in fatty acyl groups
    long_chain_pattern = Chem.MolFromSmarts("C" * 14 + "C(=O)")  # At least 14 carbons followed by C=O
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long fatty acyl chain connected"
    
    return True, "Contains sphinganine backbone with N-acyl substitution"