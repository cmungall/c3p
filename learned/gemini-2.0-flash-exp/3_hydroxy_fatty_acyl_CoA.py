"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a coenzyme A moiety linked via a thioester to the carbonyl of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified CoA pattern
    # This pattern includes the pyrophosphate, ribose, and adenine portions
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A substructure not found."

    # Thioester pattern (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, require exactly 1"

    # 3-hydroxy fatty acid pattern (-C-C(OH)-C(=O)-)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][CX4](O)[CX3](=O)")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
       return False, f"Found {len(hydroxy_matches)} 3-hydroxy groups, require exactly 1"

    # Check if the 3-hydroxy group is attached to the thioester carbonyl
    found_3hydroxy_thioester = False
    for thioester_match in thioester_matches:
        for hydroxy_match in hydroxy_matches:
          if mol.GetAtomWithIdx(thioester_match[0]).GetIdx() == mol.GetAtomWithIdx(hydroxy_match[2]).GetIdx():
            found_3hydroxy_thioester = True
            break
        if found_3hydroxy_thioester:
            break

    if not found_3hydroxy_thioester:
      return False, "Thioester not attached to 3-hydroxy carbon"

    # Check for a long carbon chain attached to the 3-hydroxy carbon and thioester carbonyl

    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds -  fatty acyl chains need to have at least 4 rotatable bonds beyond those of the core
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Fatty acid chain too short"
    
    # Minimal number of carbons should be at least 5 for a fatty acid chain plus two for the chain plus the 3-hydroxy carbon
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7 :
       return False, "Too few carbons for a 3-hydroxy fatty acyl-CoA."

    return True, "Contains CoA, thioester, 3-hydroxy fatty acid chain"