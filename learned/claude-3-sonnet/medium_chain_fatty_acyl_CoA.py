"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36034 medium-chain fatty acyl-CoA

A fatty acyl-CoA that results from the formal condensation of the thiol group of 
coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)NC(=O)CCNC(=O)CCNC(=O)S")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A substructure"
    
    # Look for acyl group
    acyl_pattern = Chem.MolFromSmarts("C(=O)")
    if not any(mol.GetSubstructMatches(acyl_pattern)):
        return False, "Missing acyl group"
    
    # Count carbon atoms in aliphatic chain
    aliphatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() < 4]
    aliphatic_chain = list(Chem.FragmentOnBonds(mol, aliphatic_atoms, set(), maxLeavingSize=20))
    chain_length = max([len(chain) for chain in aliphatic_chain])
    
    if chain_length < 6 or chain_length > 12:
        return False, "Fatty acid chain length not in medium range (6-12 carbons)"
    
    return True, "Contains Coenzyme A substructure and a medium-chain fatty acyl group"