"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl is characterized by three benzene rings connected in a 1,4 fashion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core structure identification using SMARTS
    # Represent the 1,4-diphenylbenzene core
    core_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2ccc(cc2)-c3ccccc3")
    if not mol.HasSubstructMatch(core_pattern):
      return False, "Does not contain 1,4-diphenylbenzene core"

    # 2. Ring Count Verification: Ensure exactly three benzene rings
    benzene = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene)
    if len(benzene_matches) != 3:
        return False, f"Molecule does not contain 3 benzene rings; found {len(benzene_matches)}"

    #3. Substitution check: look for common substituents such as -O, -C
    subst_pattern = Chem.MolFromSmarts("[CX4,OX2]")
    subst_matches = mol.GetSubstructMatches(subst_pattern)
    if len(subst_matches) < 3:
      return False, f"Found {len(subst_matches)} substituents, expected at least 3"
    
    #4. Molecular Weight/Size Check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbon atoms for a terphenyl"
    
    #5. Check for rotatable bonds - Terphenyls tend to be flexible
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
      return False, "Too few rotatable bonds, not flexible enough for terphenyl"
    
    #6. Check for number of attached rings.
    ring_number = rdMolDescriptors.CalcNumRings(mol)
    if ring_number < 3:
      return False, "Does not contain enough rings"

    return True, "Contains a 1,4-diphenylbenzene core with substitutions"