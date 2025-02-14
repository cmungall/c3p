"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:33569 epoxide
An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for epoxide ring pattern (O linked to 2 carbons with a single bond each)
    epoxide_pattern_1 = Chem.MolFromSmarts("[O;r3]1[C;r3][C;r3]1")
    epoxide_pattern_2 = Chem.MolFromSmarts("[O;r3]1[C;r3,H2][C;r3,H2]1")
    
    # Check if the molecule contains at least one match for the epoxide patterns
    if mol.HasSubstructMatch(epoxide_pattern_1) or mol.HasSubstructMatch(epoxide_pattern_2):
        
        # Additional checks for false positives
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        if num_rings == 0:
            return False, "No rings found in the molecule"
        
        num_3_rings = len([ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 3])
        if num_3_rings == 0:
            return False, "No 3-membered rings found in the molecule"
        
        return True, "Molecule contains an epoxide ring (a 3-membered cyclic ether)"
    else:
        return False, "No epoxide ring found in the molecule"