"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene substructure
    (benzene ring with exactly 3 chlorine substituents).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains trichlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"
    
    # For each benzene ring, check if it has exactly 3 chlorines attached
    for benzene_atoms in benzene_matches:
        chlorine_count = 0
        ring_substituents = defaultdict(int)
        
        # Check each atom in the benzene ring
        for ring_atom_idx in benzene_atoms:
            ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
            
            # Look at neighbors of ring atom
            for neighbor in ring_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine
                    chlorine_count += 1
                elif neighbor.GetIdx() not in benzene_atoms:
                    ring_substituents[neighbor.GetAtomicNum()] += 1
        
        if chlorine_count == 3:
            return True, "Found benzene ring with exactly 3 chlorine substituents"
            
    return False, "No benzene ring with exactly 3 chlorine substituents found"