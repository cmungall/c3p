"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:36689 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is a terpenoid derived from a triterpene, which has a C30 skeleton.
    The C30 skeleton may be rearranged or modified by removing skeletal atoms (usually methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for C30 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 30:
        return False, f"Molecule does not have a C30 skeleton (found {c_count} carbons)"
    
    # Look for common triterpenoid skeletons (e.g., lanostane, cycloartane, friedelane, etc.)
    skeletons = ['[C@H]1[C@@]2([C@@H](C[C@@]3([C@]1(C)CCC3(C)C)C)CCC2)C', # lanostane
                 '[C@@]12[C@@H](C[C@@H](C1)C)CCC3=[C]2CCC4[C@]3(CCC(=C4)C)C', # cycloartane
                 '[C@@]12[C@@H](C[C@H](C1)C)CCC3=[C]2[C@H](CC4[C@@]3(CCC(=C4)C)C)C', # friedelane
                 # Add more skeletons as needed
                 ]
    for skeleton in skeletons:
        skeleton_mol = Chem.MolFromSmarts(skeleton)
        if mol.HasSubstructMatch(skeleton_mol):
            return True, f"Contains {skeleton_mol.GetProp('_Name')} skeleton, a common triterpenoid scaffold"
    
    # Check for rearranged/modified C30 skeleton
    rearranged_pattern = Chem.MolFromSmarts('[C@H]1[C@]2([C@H](C[C@@]3([C@]1(C)CCC3(C)C)C)CCC2)C' # lanostane with modified/removed methyl groups
                                             '|[C@@]12[C@@H](C[C@@H](C1)C)CCC3=[C]2CCC4[C@]3(CCC(=C4)C)C' # cycloartane with modified/removed methyl groups
                                             '|[C@@]12[C@@H](C[C@H](C1)C)CCC3=[C]2[C@H](CC4[C@@]3(CCC(=C4)C)C)C') # friedelane with modified/removed methyl groups
    if mol.HasSubstructMatch(rearranged_pattern):
        return True, "Contains a rearranged/modified C30 skeleton derived from a common triterpenoid scaffold"
    
    return False, "No evidence of a triterpenoid scaffold"