"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are derived from prostanoic acid, a C20 fatty acid with a cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for cyclopentane ring with attached chains
    cyclopentane_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1")  #Basic 5 member ring
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    #2. Check for a carboxyl acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    #3. Check for 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Not 20 carbons, found {c_count}"
    
    #4. Check for at least one double bond
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond found"

    # 5. Check for 4-6 oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4 or o_count > 6:
        return False, f"Oxygen count is not between 4 and 6, found {o_count}"

    return True, "Meets basic prostaglandin criteria"