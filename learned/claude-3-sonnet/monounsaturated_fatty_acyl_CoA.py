"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:49823 monounsaturated fatty acyl-CoA

A monounsaturated fatty acyl-CoA is any unsaturated fatty acyl-CoA in which the fatty acyl chain
contains one carbon-carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"
    
    # Look for fatty acyl chain (long carbon chain)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "No fatty acyl chain found"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acyl chain too short"
    
    # Look for exactly one double bond in the chain
    unsaturation_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    unsaturation_matches = mol.GetSubstructMatches(unsaturation_pattern)
    if len(unsaturation_matches) != 1:
        return False, f"Found {len(unsaturation_matches)} double bonds, expected 1"
    
    return True, "Contains a CoA group attached to a monounsaturated fatty acyl chain"