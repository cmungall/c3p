"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
Definition: An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define the coenzyme A (CoA) substructure pattern
    # Simplified pattern for the CoA moiety
    coa_smarts = """
    N[C@@H](CCNC(=O)[C@@H](O)C(C)(C)COP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@@H](O)[C@H]1OP(=O)(O)O)C(=O)O
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not coa_pattern:
        return False, "Invalid CoA SMARTS pattern"

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    # Define the thioester linkage pattern
    # Thioester linkage between fatty acyl chain and CoA: C(=O)-S-C
    thioester_smarts = "C(=O)SCCNC(=O)"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not thioester_pattern:
        return False, "Invalid thioester SMARTS pattern"

    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage to CoA not found"

    # Define the 3-oxo-fatty acyl chain pattern
    # Fatty acyl chain with a keto group at the 3-position
    # Pattern: O=C-C-C(=O)-C
    three_oxo_acyl_smarts = "C(=O)CC(=O)C"
    three_oxo_acyl_pattern = Chem.MolFromSmarts(three_oxo_acyl_smarts)
    if not three_oxo_acyl_pattern:
        return False, "Invalid 3-oxo-fatty acyl SMARTS pattern"

    # Check for 3-oxo-fatty acyl chain attached via thioester linkage
    if not mol.HasSubstructMatch(three_oxo_acyl_pattern):
        return False, "3-oxo-fatty acyl chain not found"

    # Optionally, check the length of the fatty acyl chain
    # Find the match for the acyl chain
    acyl_matches = mol.GetSubstructMatches(three_oxo_acyl_pattern)
    if not acyl_matches:
        return False, "3-oxo-fatty acyl chain not found"

    for match in acyl_matches:
        # Get the acyl chain atoms from the match
        acyl_chain_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Count the number of carbons in the acyl chain
        c_count = sum(1 for atom in acyl_chain_atoms if atom.GetAtomicNum() == 6)
        if c_count < 5:  # Minimum length to be considered a fatty acid
            return False, "Fatty acyl chain too short"

    return True, "Molecule is a 3-oxo-fatty acyl-CoA"