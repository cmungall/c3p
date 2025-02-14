"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is characterized by a glycerol backbone with three esterified fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern allowing for stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 3 ester groups (-C(=O)O-)
    # We want exactly three ester bonds, consistent with triglyceride structure
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, expected 3"

    # Recognize fatty acid chains (saturated/unsaturated)
    # Allow for chains with extensions reflecting unsaturation (alkene, alkyne presence)
    long_chain_min_length = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0]")
    fatty_acid_matches = mol.GetSubstructMatches(long_chain_min_length)
    
    # Check for enough potential fatty chains
    if len(fatty_acid_matches) < 3:
        return False, f"Missing sufficient fatty acid-like chains, detected {len(fatty_acid_matches)}"

    # Broader molecular weight range reflecting natural triglyceride variability
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300:
        return False, "Molecular weight is quite low, typically >300 Da for triglycerides"

    return True, "SMILES represents a triglyceride with a glycerol backbone and three ester groups"