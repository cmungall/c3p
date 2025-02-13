"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: CHEBI:15948 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid occurring in bile,
    usually with an amide linkage to glycine or taurine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cholanic acid backbone (tetracyclic steroid)
    cholanic_pattern = Chem.MolFromSmarts("[C@@H]3[C@@]2([C@@H]([C@]1([C@H](CC[C@]4([C@@H](CC(=O)O)C[C@H]4CC4)CC1)C)CC[C@@H]2C[C@@H]3O)C")
    if not mol.HasSubstructMatch(cholanic_pattern):
        return False, "No cholanic acid backbone found"

    # Check for 5beta configuration (cis junction rings B/C)
    stereo_config = Chem.FindMolChiralUnspecifiedUnknown(mol, includeUnknown=True)
    if any(code == Chem.ChiralType.CHI_TETRAHEDRAL_CIS for code in stereo_config):
        return False, "Not 5beta configuration (trans B/C ring junction)"

    # Look for hydroxy groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxy groups found"

    # Look for carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Look for amide linkage to glycine or taurine (optional)
    amide_pattern = Chem.MolFromSmarts("[NX3H2]-[CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        return True, "Contains cholanic acid backbone with hydroxy groups and amide linkage"

    return True, "Contains cholanic acid backbone with hydroxy groups"