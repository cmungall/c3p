"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:27092 porphyrins
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_porphyrin(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins are natural pigments containing a fundamental skeleton of four pyrrole nuclei
    united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for porphyrin macrocycle
    porphyrin_pattern = Chem.MolFromSmarts("c1nc2nc3nc4nc1ccc4cc3cc2")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin macrocycle found"
    
    # Check for 4 pyrrole rings
    pyrrole_pattern = Chem.MolFromSmarts("c1cc[nH]c1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) != 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole rings, need exactly 4"
    
    # Check for methine groups connecting pyrroles
    methine_pattern = Chem.MolFromSmarts("C=C")
    methine_matches = mol.GetSubstructMatches(methine_pattern)
    if len(methine_matches) != 4:
        return False, f"Found {len(methine_matches)} methine groups, need exactly 4"
    
    # Check for aromaticity
    if not mol.GetAromaticRingInfo():
        return False, "Porphyrin macrocycle is not aromatic"
    
    # Check for substituents (metals, side chains, etc.)
    # This is not strictly required, but porphyrins usually have substituents
    substituents = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [5, 6, 7, 8]]
    if not substituents:
        return True, "Contains bare porphyrin macrocycle"
    else:
        return True, f"Contains porphyrin macrocycle with substituents: {', '.join(set(substituents))}"