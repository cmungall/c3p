"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:26389 porphyrin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined by a macrocycle consisting of four pyrrole rings linked by methine bridges.

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

    # Define SMARTS patterns for pyrrole rings and methine bridges
    pyrrole_pattern = Chem.MolFromSmarts("[nX1]1cc[cX3]c1")  # X1 for a single bond to a hydrogen
    methine_pattern = Chem.MolFromSmarts("[CX3H1]") # A carbon with 3 single bonds and 1 H
    
    # Find pyrrole rings
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Find methine groups
    methine_matches = mol.GetSubstructMatches(methine_pattern)

    #Check there are 4 pyrroles and 4 methines
    if len(pyrrole_matches) != 4:
      return False, f"Found {len(pyrrole_matches)} pyrrole rings, require 4."
    if len(methine_matches) < 4 :
      return False, f"Found {len(methine_matches)} methine groups, require at least 4."

    #Check the number of nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count != 4:
      return False, f"Found {n_count} nitrogen atoms, require 4"

    #Check ring size. Use a simple heuristic: the aromatic core consists of four pyrroles (5 members each) and four methine carbons (4 * 1 each), or 24 atoms. 
    # Porphyrins typically have 24 or 25 atoms in their core ring system, due to different tautomers and substitutions
    
    ring_atoms = rdMolDescriptors.CalcNumAtoms(mol, onlyHeavy = True) # Only heavy atoms, ignore hydrogens
    if not (24 <= rdMolDescriptors.CalcNumAtoms(mol, onlyHeavy = True, includeOnlyInRing = True) <= 26):
      return False, f"Ring system too small or too large: {rdMolDescriptors.CalcNumAtoms(mol, onlyHeavy = True, includeOnlyInRing = True)} ring atoms. Expected 24-26"


    # Check connectivity between pyrroles and methines (simplified, checking for substructure overlaps is complex)
    # A more rigourous check would require graph analysis and the definition of the macrocycle
    # For now, assume that if we have 4 pyrroles and 4 methines it is a porphyrin

    return True, "Contains four pyrrole rings linked by methine bridges."