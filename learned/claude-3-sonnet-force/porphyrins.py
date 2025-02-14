"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:27026 porphyrins
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_porphyrin(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is a macrocyclic compound containing four pyrrole rings linked
    by methine bridges.

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
    
    # Look for porphyrin macrocycle pattern
    porphyrin_pattern = Chem.MolFromSmarts("[N,n]1c2ccc[nH]c2c3ccc[nH]c3c4ccc[nH]c4c1")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin macrocycle found"
    
    # Count pyrrole rings
    pyrrole_pattern = Chem.MolFromSmarts("c1ccc[nH]1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    num_pyrroles = len(pyrrole_matches)
    if num_pyrroles != 4:
        return False, f"Found {num_pyrroles} pyrrole rings, need exactly 4"
    
    # Check for methine bridges
    methine_pattern = Chem.MolFromSmarts("C(=[N,n])")
    methine_matches = mol.GetSubstructMatches(methine_pattern)
    num_methines = len(methine_matches)
    if num_methines != 4:
        return False, f"Found {num_methines} methine bridges, need exactly 4"
    
    # Check for planar macrocycle
    planar = AllChem.IsPlanar(mol)
    if not planar:
        return False, "Porphyrin macrocycle is not planar"
    
    # Check for metal complexation
    metal_pattern = Chem.MolFromSmarts("[Mg,Fe,Co,Ni,Cu,Zn,Pd,Pt]")
    metal_matches = mol.GetSubstructMatches(metal_pattern)
    if len(metal_matches) > 0:
        return True, "Contains porphyrin macrocycle with metal complexation"
    else:
        return True, "Contains porphyrin macrocycle without metal complexation"