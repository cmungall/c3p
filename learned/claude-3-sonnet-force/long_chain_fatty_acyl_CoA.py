"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35621 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA backbone pattern
    CoA_pattern = Chem.MolFromSmarts("C1OC(COP(=O)([O-])OP(=O)([O-])OCC2OC(N3C=NC4=C3N=CN=C4N)C(O)C2O)C(O)C1"
    if not mol.HasSubstructMatch(CoA_pattern):
        return False, "No CoA backbone found"

    # Look for fatty acyl chain pattern (long carbon chain with terminal carbonyl)
    fatty_acyl_pattern = Chem.MolFromSmarts("CCC(=O)CCCCCCCCCCCC")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    if len(fatty_acyl_matches) != 1:
        return False, f"Found {len(fatty_acyl_matches)} fatty acyl chains, need exactly 1"

    # Count carbons in the fatty acyl chain
    fatty_acyl_chain = Chem.DeleteSubstructs(mol, CoA_pattern)
    n_carbons = sum(1 for atom in fatty_acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 13 or n_carbons > 22:
        return False, f"Found {n_carbons} carbons in fatty acyl chain, should be between 13 and 22"

    # Count rotatable bonds in the fatty acyl chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(fatty_acyl_chain)
    if n_rotatable < 10:
        return False, "Fatty acyl chain too short or rigid"

    return True, "Contains CoA backbone with a long-chain fatty acyl group"