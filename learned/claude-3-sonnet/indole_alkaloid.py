"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:38256 indole alkaloid
An alkaloid containing an indole skeleton.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for indole substructure
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cnc2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole substructure found"
    
    # Check for nitrogen atom (alkaloid criteria)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7) == 0:
        return False, "No nitrogen atom found (alkaloid criteria)"
    
    # Check molecular weight (typical for alkaloids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, "Molecular weight too high for alkaloid"
    
    # Check for aromatic rings (typical for indole alkaloids)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, "Too few aromatic rings for indole alkaloid"
    
    # Additional checks or structural patterns can be added here
    
    return True, "Contains an indole substructure and meets alkaloid criteria"