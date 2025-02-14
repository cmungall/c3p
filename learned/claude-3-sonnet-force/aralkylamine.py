"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: CHEBI:32935 aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for aromatic ring(s)
    aromatic_rings = mol.GetAromaticRings()
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    # Look for primary/secondary amine group(s)
    amine_pattern = Chem.MolFromSmarts("[NH2,NH1]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No primary or secondary amine groups found"
    
    # Check if amine group is directly attached to aromatic ring
    for amine_idx in amine_matches:
        amine_atom = mol.GetAtomWithIdx(amine_idx)
        for neighbor in amine_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Contains an alkylamine with aromatic substituent"
    
    # Amine group(s) not directly attached to aromatic ring
    return False, "Amine group(s) not directly attached to aromatic ring"