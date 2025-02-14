"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:59382 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4_(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a fatty acyl-CoA with a long carbon chain arising from deprotonation
    of the phosphate and diphosphate OH groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("C(C(C(=O)NCCC(=O)NCCS)O)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for fatty acid chain (long carbon chain with carboxyl group)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CCCCCCC")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, f"Found {len(fatty_acid_matches)} fatty acid chains, need exactly 1"
    
    # Check for deprotonated phosphate/diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])([O-])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} deprotonated phosphate groups, need exactly 2"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acid chain too short"
    
    return True, "Contains CoA backbone with a long fatty acid chain and deprotonated phosphate/diphosphate groups"