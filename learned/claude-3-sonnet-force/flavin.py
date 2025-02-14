"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: CHEBI:28363 flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavins are derivatives of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
    skeleton, with a substituent on the 10 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoalloxazine backbone pattern with dimethyl groups on the benzene ring
    backbone_pattern = Chem.MolFromSmarts("c1nc2c(nc3c(=O)[nH]c(=O)nc23)c(C)c(C)c1")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Missing dimethylisoalloxazine backbone"
    
    # Look for substituent at the 10 position
    # Use a more general pattern to match a wide range of substituents
    substituent_pattern = Chem.MolFromSmarts("[C,N][N]1C=2C(=NC3=C1C=CC=C3)C(NC(N2)=O)=O[*]")
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "No substituent at the 10 position"
    
    # Check for atoms and rings characteristic of flavins
    n_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if 6 not in n_atoms or 7 not in n_atoms or 8 not in n_atoms:
        return False, "Missing required atoms (C, N, O)"
    
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(x) for x in ring_info.AtomRings()]
    if 5 not in ring_sizes or 6 not in ring_sizes:
        return False, "Missing required 5-membered and 6-membered rings"
    
    return True, "Contains dimethylisoalloxazine backbone with substituent on the 10 position"