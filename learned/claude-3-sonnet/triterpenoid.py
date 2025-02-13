"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:36691 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is a terpenoid derived from a triterpene with a C30 skeleton,
    which may be rearranged or modified by removal of skeletal atoms (e.g. methyl groups).

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
    
    # Check for 30 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 30:
        return False, f"Found {c_count} carbon atoms, expected 30 for triterpenoid"
    
    # Check for tetracyclic structure
    cycle_list = mol.GetRingInfo().AtomRings()
    if len(cycle_list) < 4:
        return False, "Less than 4 rings, expected tetracyclic structure"
    
    # Check for long carbon chains (lipophilic moieties)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long carbon chains found, expected lipophilic moieties"
    
    # Check for triterpene skeleton
    triterpene_pattern = Chem.MolFromSmarts("[C@@H]2[C@@H]1[C@H]([C@@H]([C@H](C1)[C@@H]2C)C)C")
    if not mol.HasSubstructMatch(triterpene_pattern):
        return False, "No triterpene skeleton found"
    
    # Passed all checks, classify as triterpenoid
    return True, "Contains a tetracyclic structure with 30 carbon atoms and lipophilic moieties, derived from a triterpene"