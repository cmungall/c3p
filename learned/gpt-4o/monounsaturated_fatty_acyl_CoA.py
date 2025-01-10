"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    Such a molecule contains a Coenzyme A structure and a fatty acyl chain with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for extended Coenzyme A structure with charge considerations
    coa_pattern = Chem.MolFromSmarts(
        "[NX3][CX3](=O)[CNH][CX3](=O)[C@H](O)[C](C)(C)[COP](=O)([O])O[PX4](=O)O[CX2H][C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(nc[nH]c12)"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found or incorrect charge state"

    # Check for exactly one carbon-carbon double bond in the acyl chain
    db_pattern = Chem.MolFromSmarts("C=C")
    db_matches = mol.GetSubstructMatches(db_pattern)
    
    # Count double bonds that are part of fatty acyl chain, ensuring they are not within the CoA structure
    fatty_acyl_db_counter = 0
    for match in db_matches:
        if any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6 for atom_idx in match):
            # Ensure the double bond is not within the commonly known structure of CoA (using index heuristics)
            if all(mol.GetAtomWithIdx(atom_idx).GetIdx() < 40 for atom_idx in match):
                fatty_acyl_db_counter += 1

    # Verify there is exactly one double bond in the acyl chain
    if fatty_acyl_db_counter != 1:
        return False, f"Expected 1 carbon-carbon double bond, found {fatty_acyl_db_counter}"

    return True, "Contains Coenzyme A structure and a fatty acyl chain with exactly one carbon-carbon double bond"