"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:72737 trans-2-enoyl-CoA

A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA that results from the formal condensation
of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A backbone
    coenzyme_A_pattern = Chem.MolFromSmarts("C(C)(C)(CO)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coenzyme_A_pattern):
        return False, "No coenzyme A backbone found"

    # Look for trans-2-enoyl group
    trans_2_enoyl_pattern = Chem.MolFromSmarts("[CX3](/C=C/[CX3]([CX3])=[OX1])")
    if not mol.HasSubstructMatch(trans_2_enoyl_pattern):
        return False, "No trans-2-enoyl group found"

    # Look for fatty acid chain (long carbon chain)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)} matches"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for trans-2-enoyl-CoA"

    # Count carbons, oxygens, nitrogens, and sulfurs
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    if c_count < 20:
        return False, "Too few carbons for trans-2-enoyl-CoA"
    if o_count != 10:
        return False, "Must have exactly 10 oxygens"
    if n_count != 3:
        return False, "Must have exactly 3 nitrogens"
    if s_count != 1:
        return False, "Must have exactly 1 sulfur"

    return True, "Contains coenzyme A backbone with trans-2-enoyl group and fatty acid chain"