"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for unsaturation (C=C or C#C bonds)
    unsaturated = any(bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
                      for bond in mol.GetBonds())
    if not unsaturated:
        return False, "No C=C or C#C bonds found (saturated)"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not chain_matches:
        return False, "Carbon chain too short"
    
    # Passed all checks, classify as unsaturated fatty acid
    return True, "Contains at least one C=C or C#C bond and a carboxylic acid group"