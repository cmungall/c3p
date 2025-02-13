"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:53236 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is an oxo fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.

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

    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])O)n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Look for 3-oxo group on fatty acid chain
    oxo_pattern = Chem.MolFromSmarts("CC(=O)CC")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group on fatty acid chain"

    # Look for long carbon chain (>= 6 carbons)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "Fatty acid chain too short"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Fatty acid chain too short"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for 3-oxo-fatty acyl-CoA"
    if o_count < 10:
        return False, "Too few oxygens for 3-oxo-fatty acyl-CoA"

    return True, "Contains CoA substructure and 3-oxo group on fatty acid chain"