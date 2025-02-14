"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36345 very long-chain fatty acyl-CoA
A fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(CO[P@](=O)(O)O[P@](=O)(O)OC[C@H]1[C@@H]([C@@H](O)[C@H](O)OP(O)(O)=O)O[C@@H]1n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCSC(=O)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA backbone"
    
    # Look for long fatty acid chain (>22 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]([CX3])[CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Fatty acid chain too short (<=C22)"
    
    # Count double bonds in fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    n_double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Check for any other functional groups
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6, 8, 7, 15, 16)) > 0:
        return False, "Contains other functional groups besides C, H, O, N, P, S"
    
    # Classify based on chain length and double bonds
    chain_length = rdMolDescriptors.CalcNumAliphaticCarbons(mol)
    reason = f"Very long-chain fatty acyl-CoA with {chain_length} carbons and {n_double_bonds} double bonds"
    
    return True, reason