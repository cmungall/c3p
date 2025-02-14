"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:83250 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA results from the formal condensation of the thiol group of coenzyme A 
    with the carboxy group of any 3-oxo-fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for coenzyme A substructure
    coa_pattern = Chem.MolFromSmarts("[N]C(CCNC(=O)CC[N]C([C@@H](O)[C@H](C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)=O)n1cnc2c(N)ncnc12))(O)=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A substructure found"
    
    # Look for oxo group at position 3
    oxo_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No 3-oxo group found"
    
    # Check for fatty acid chain (long carbon chain)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Fatty acid chain too short"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Fatty acid chain too short"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for a fatty acid"
    if o_count < 4:
        return False, "Too few oxygens for 3-oxo-fatty acyl-CoA"
    
    return True, "Contains coenzyme A substructure and a 3-oxo group with a fatty acid chain"