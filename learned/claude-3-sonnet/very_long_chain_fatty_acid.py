"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:36979 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Count carbon chain length
    carbon_chain = Chem.MolFromSmarts("[CH4]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]")
    chain_matches = mol.GetSubstructMatches(carbon_chain)
    chain_lengths = [len(match) for match in chain_matches]
    max_chain_length = max(chain_lengths, default=0)

    if max_chain_length <= 22:
        return False, f"Longest carbon chain has {max_chain_length} atoms, must be > 22"

    # Check for unsaturations
    unsaturated_chain = Chem.MolFromSmarts("[CH2]=[CH]~[CH2]~[CH2]")
    unsaturated_matches = mol.GetSubstructMatches(unsaturated_chain)
    unsaturated_lengths = [len(match) for match in unsaturated_matches]
    max_unsaturated_length = max(unsaturated_lengths, default=0)

    if max_unsaturated_length > 4:
        return False, "Too many unsaturations for fatty acid"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if "Br" in formula or "Cl" in formula or "F" in formula or "I" in formula:
        return False, "Halogenated compounds not allowed"

    return True, "Carbon chain length > 22, with carboxylic acid group"