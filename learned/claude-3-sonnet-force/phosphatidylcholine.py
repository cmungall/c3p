"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: CHEBI:17556 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine with two acyl substituents at positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphocholine group (-OP(=O)([O-])OCC[N+](C)(C)C)
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
    
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for long carbon chains (fatty acid chains) attached to esters
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = []
    for ester_idx in ester_matches:
        ester_atom = mol.GetAtomWithIdx(ester_idx)
        for neighbor in ester_atom.GetNeighbors():
            chain_match = mol.GetSubstructMatches(fatty_acid_pattern, atomIdx=neighbor.GetIdx())
            fatty_acid_matches.extend(chain_match)
    if len(set(fatty_acid_matches)) < 2:
        return False, "Missing fatty acid chains"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylcholines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylcholine"

    return True, "Contains glycerol backbone with 2 fatty acid chains and phosphocholine group"

# Example usage
smiles = "CCCCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC/C=C\C/C=C\C/C=C\CC"
is_pc, reason = is_phosphatidylcholine(smiles)
print(f"Is phosphatidylcholine? {is_pc}")
print(f"Reason: {reason}")