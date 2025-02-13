"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is defined as any 3-hydroxy steroid whose skeleton is closely related
    to cholestan-3-ol, possibly with additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a 3-hydroxy steroid core: steroid skeleton with a hydroxyl group at C3
    steroid_core_pattern = Chem.MolFromSmarts("C1CC[C@H]2[C@@H]3CC[C@]4(C)CC[C@H](O)C4[C@]3(C)CC[C@]12C")
    
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain a 3-hydroxy steroid core structure"
        
    # Additional check for optional side chain
    # This will be highly context-specific, so the minimum steroid pattern check is typically sufficient for sterols
    
    return True, "Contains a 3-hydroxy steroid structure related to cholestan-3-ol"