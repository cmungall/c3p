"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35571 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA results from the condensation of coenzyme A with a medium-chain fatty acid.

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

    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H](OP(O)(O)=O)O1)n1cnc2c(N)ncnc12)OP(OP(OCC(C)(C(=O)NCCC(=O)NCCS)O)(=O)O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Look for medium-chain fatty acid (6-12 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No medium-chain fatty acid found"
    
    # Count carbons in fatty acid chain
    fatty_acid_carbons = [len(Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, match, allHsExplicit=False, isomericSmiles=True)).GetAtoms(onlyExplicit=True)) for match in fatty_acid_matches]
    if not any(6 <= n <= 12 for n in fatty_acid_carbons):
        return False, "Fatty acid chain length not in medium range (6-12 carbons)"

    # Check for ester bond
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Must have exactly one ester bond"

    return True, "Contains CoA backbone with medium-chain fatty acid attached via ester bond"