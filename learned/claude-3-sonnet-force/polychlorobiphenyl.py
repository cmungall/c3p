"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: CHEBI:38138 polychlorobiphenyl

A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for biphenyl core pattern
    biphenyl_pattern = Chem.MolFromSmarts("c1ccc(-c2ccccc2)cc1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"
    
    # Look for chlorine atoms attached to the biphenyl core
    cl_pattern = Chem.MolFromSmarts("c1ccc(-c2ccc(Cl)cc2)cc1")
    cl_matches = mol.GetSubstructMatches(cl_pattern)
    n_cl_atoms = len(cl_matches)
    if n_cl_atoms < 2 or n_cl_atoms > 10:
        return False, f"Found {n_cl_atoms} chlorine atoms, must be between 2 and 10"
    
    # Check molecular weight range (200-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside expected range"
    
    # Check for other heteroatoms (should only have C, H, and Cl)
    atoms = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    allowed_atoms = set([1, 6, 17])  # H, C, Cl
    if not atoms.issubset(allowed_atoms):
        return False, "Found heteroatoms other than C, H, and Cl"
    
    return True, "Contains a biphenyl core with 2-10 chlorine atoms attached"