"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:36654 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene (C15 skeleton),
    which may have rearranged or modified by the removal of one or more methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Found {c_count} carbon atoms, sesquiterpenoids must have exactly 15"
    
    # Look for characteristic sesquiterpene skeleton patterns
    # This is a non-exhaustive list, more patterns could be added
    patterns = ['CC(C)CCCC(C)C', 'CC(C)CCCCCC(C)C', 'CC(C)CC(C)CC(C)', 'CC(C)CCC(C)(C)C', 'CC(C)CCCC(C)(C)C']
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return True, "Found sesquiterpene skeleton pattern"
    
    # If no skeleton pattern matched, check other criteria
    
    # Count rings - terpenoids typically have 1-4 rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 1 or n_rings > 4:
        return False, f"Found {n_rings} rings, sesquiterpenoids typically have 1-4 rings"
    
    # Count rotatable bonds - terpenoids have few rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for a sesquiterpenoid"
    
    # If all criteria matched, assume it is a sesquiterpenoid
    return True, "Meets criteria for sesquiterpenoid (15 carbons, ring/rotatable bond counts, potential skeleton pattern)"