"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    Azoles are monocyclic 5-membered heterocycles containing at least one nitrogen and
    possibly other heteroatoms (N, S, O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find 5-membered rings with at least 1 Nitrogen
    # This SMARTS pattern matches any 5-membered ring with at least one N and optionally other heteroatoms (O, S, N)
    azole_pattern = Chem.MolFromSmarts("[NX2,NX3;R5][CX3,CX2,NX2,OX2,SX2;R5]~[CX3,CX2,NX2,OX2,SX2;R5]~[CX3,CX2,NX2,OX2,SX2;R5]~[CX3,CX2,NX2,OX2,SX2;R5]") # 5-membered ring with at least one N
    
    ring_matches = mol.GetSubstructMatches(azole_pattern)

    if not ring_matches:
          return False, "No azole ring found"

    #Check that all atoms in the ring have ring bond count = 2
    for match in ring_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetDegree() > 3: # If total degree >3 it is part of fused ring. 2 ring bonds max + 1 substituent is ok
              return False, "Found fused ring or too many bonds to ring atom"
            if atom.IsInRing() and atom.GetRingBondCount() != 2: # check for monocyclic ring
              return False, "ring atom not part of monocyclic ring"

    
    return True, "Contains a monocyclic 5-membered heterocycle with at least one nitrogen."