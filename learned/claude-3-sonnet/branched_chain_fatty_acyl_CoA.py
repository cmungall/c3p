"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36565 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA that results from the formal
    condensation of the thiol group of coenzyme A with the carboxy group of any branched-chain
    fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for branched fatty acid chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4]([CX4,CX3])([CX4,CX3])")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No branched fatty acid chain found"
    
    # Check for ester bond between CoA and fatty acid
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"
    
    return True, "Contains CoA backbone and branched fatty acid chain joined by ester bond"