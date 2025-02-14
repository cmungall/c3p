"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine consists of a sphinganine backbone with a fatty acyl group attached to the nitrogen.

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

    # Revised Sphinganine backbone pattern, does not define the chain length.
    sphinganine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CHX4]([CH2X4]O)N")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"

    # N-acyl bond pattern (N-C=O), where the C=O is attached to a long chain
    acyl_pattern = Chem.MolFromSmarts("N[CX3](=[OX1])[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No N-acyl group found"
    
    # Verify that there is only one such match (one N-acyl group).
    if len(acyl_matches) != 1:
         return False, "More than one N-acyl group"

    # Check for fatty acid chain on the acyl group
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
      return False, f"Missing fatty acyl chain"

    # Check for a minimum of rotatable bonds for the fatty chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 7:
        return False, "Chains too short"

    # Count Carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14: # a minimum number of carbons in sphinganine backbone
        return False, "Too few carbons for a N-acylsphinganine"
    
    # Count Oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
      return False, "Too few oxygens"


    return True, "Contains sphinganine backbone with a fatty acyl group attached to the nitrogen"