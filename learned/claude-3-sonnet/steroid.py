"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:17240 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as 'Any of naturally occurring compounds and synthetic analogues, based on the cyclopenta[a]phenanthrene carbon skeleton, partially or completely hydrogenated; there are usually methyl groups at C-10 and C-13, and often an alkyl group at C-17. By extension, one or more bond scissions, ring expansions and/or ring contractions of the skeleton may have occurred. Natural steroids are derived biogenetically from squalene which is a triterpene.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopenta[a]phenanthrene skeleton
    steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@H]3CC[C@@]4([C@H](CC[C@]4([C@H]3CC[C@@]21C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No cyclopenta[a]phenanthrene skeleton found"
    
    # Check for methyl groups at C-10 and C-13
    c10_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@H]3CC[C@@]4([C@H](CC[C@]4([C@H]3CC[C@]21C)C)[C@H](C)C)C)C")
    c13_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@H]3CC[C@@]4([C@H](CC[C@]4([C@H]3CC[C@]21C)C)C)[C@@H](C)C)C")
    if not (mol.HasSubstructMatch(c10_pattern) and mol.HasSubstructMatch(c13_pattern)):
        return False, "Missing methyl groups at C-10 and/or C-13"
    
    # Check for alkyl group at C-17 (optional)
    c17_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@]2([C@H]3CC[C@@]4([C@H](CC[C@]4([C@H]3CC[C@]21C)C)CC)C)C")
    has_c17_alkyl = mol.HasSubstructMatch(c17_pattern)
    
    # Check for bond scissions, ring expansions/contractions
    # Use Bemis-Murcko scaffold to identify skeleton modifications
    scaffold = AllChem.MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C5CCC6C7CCC8C9CCCC%91%92%93%94%95%96%97%98%99")
    if not scaffold.HasSubstructMatch(scaffold_pattern):
        return False, "Skeleton modified by bond scissions, ring expansions/contractions"
    
    # Check molecular weight - steroids typically 200-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical steroid range"
    
    reason = "Contains cyclopenta[a]phenanthrene skeleton with methyl groups at C-10 and C-13"
    if has_c17_alkyl:
        reason += ", and an alkyl group at C-17"
    
    return True, reason