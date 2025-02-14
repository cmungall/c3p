"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core isoflavone structure using SMARTS
    core_smarts = "c1cc2c(=O)c(-c3ccccc3)oc2cc1"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core isoflavone skeleton not found"

    # Check that there are no additional fused rings
    fused_ring_pattern = Chem.MolFromSmarts("c1cc2c3ccccc3oc2cc1")
    fused_matches = mol.GetSubstructMatches(fused_ring_pattern)
    if len(fused_matches) != 1:
      return False, "Additional fused rings found"

    # Check that the carbonyl group is present
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
      return False, "Core carbonyl not found"
    
    # Check that there is exactly 1 aryl group in position 3
    aryl_pattern = Chem.MolFromSmarts("c1(-c2ccccc2)ccccc1")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern)
    if len(aryl_matches) !=1:
       return False, "Exactly one aryl group in position 3 is required"

    # Further checks (optional, might be too strict for certain derivatives)
    # Count number of carbons, oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 15:
        return False, "Too few carbons for isoflavone"
    if o_count < 2:
        return False, "Too few oxygens for isoflavone"


    return True, "Matches isoflavone skeleton"