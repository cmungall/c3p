"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:33309 fatty acyl-CoA
An acyl-CoA that results from the formal condensation of the thiol group of 
coenzyme A with the carboxy group of any fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A backbone pattern
    coa_pattern = Chem.MolFromSmarts("[C@H]1[N@+]2=CN=C(N)N=C2N1C3OC(COP(=O)([O-])OP(=O)([O-])OP(=O)(O)O)C(OP(=O)([O-])O)C3O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found"

    # Look for thioester (-C(=O)-S-) linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"

    # Check for fatty acid chain (long carbon chain attached to thioester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chain too short to be fatty acid"

    return True, "Contains coenzyme A backbone with a fatty acid chain attached via a thioester bond"