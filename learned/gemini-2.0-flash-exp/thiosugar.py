"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"
    
     # Check for at least one oxygen atom, or for sulfur replacing an oxygen in the ring.
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygen_atoms and not any(atom.IsInRing() for atom in sulfur_atoms):
        return False, "No oxygen atoms found, and sulfur not in ring (not a carbohydrate derivative)"
    
    # Detect ring systems
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    if not rings:
      # If no rings found, check for -SR type structures
      for s_atom in sulfur_atoms:
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6: # carbon
                # Ensure the S is not part of sulfone/sulfonamide or similar by checking neighbors of the carbon
                is_sulfone_or_similar = False
                for carbon_neighbor in neighbor.GetNeighbors():
                  if carbon_neighbor.GetAtomicNum() == 16 and carbon_neighbor != s_atom: # check if attached to another S
                      is_sulfone_or_similar = True
                      break
                  elif carbon_neighbor.GetAtomicNum() == 8: # oxygen, also check if attached to =O
                      for oxy_neighbor in carbon_neighbor.GetNeighbors():
                          if oxy_neighbor.GetAtomicNum() == 16 and oxy_neighbor == s_atom:
                             is_sulfone_or_similar = True
                             break
                      if is_sulfone_or_similar:
                          break
                if not is_sulfone_or_similar:
                   return True, "Sulfur attached to a carbon, potentially -SR"
      return False, "No rings found, and sulfur not attached to C in a suitable way"

    #Check for rings of size 5 or 6
    valid_rings = [ring for ring in rings if len(ring) == 5 or len(ring) == 6]
    if not valid_rings:
      return False, "No 5- or 6-membered rings found"
    
    for ring in valid_rings:
         # Check if the ring contains mostly carbons, and one O or one S
         carbon_count = 0
         heteroatom_count = 0
         has_sulfur = False
         has_oxygen = False
         for atom_idx in ring:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetAtomicNum() == 6:
                 carbon_count += 1
             elif atom.GetAtomicNum() == 8:
                 heteroatom_count += 1
                 has_oxygen = True
             elif atom.GetAtomicNum() == 16:
                 heteroatom_count += 1
                 has_sulfur = True
            
         if carbon_count + heteroatom_count != len(ring) or heteroatom_count > 1:
            continue # Not carbohydrate-like
         
         # Check if an oxygen in the ring is replaced by sulfur or if any hydroxyl group is replaced by sulfur
         if has_sulfur:
           if (has_oxygen or carbon_count > 0):
             return True, "Sulfur in ring replacing oxygen"
             
           
         # Check for -SR substitutions for -OH groups outside the ring, and that the sulfur does not belong to a sulfonamide/sulfone
         for atom_idx in ring:
              atom = mol.GetAtomWithIdx(atom_idx)
              if atom.GetAtomicNum() == 6:
                   for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0: # check OH- group
                         for s_atom in sulfur_atoms:
                           if neighbor in s_atom.GetNeighbors() or any(n.GetAtomicNum() == 6 for n in s_atom.GetNeighbors()):
                              # Ensure the S is not part of sulfone/sulfonamide or similar
                              is_sulfone_or_similar = False
                              for carbon_neighbor in atom.GetNeighbors():
                                  if carbon_neighbor.GetAtomicNum() == 16 and carbon_neighbor != s_atom: #check if attached to another S
                                      is_sulfone_or_similar = True
                                      break
                                  elif carbon_neighbor.GetAtomicNum() == 8: # oxygen, also check if attached to =O
                                      for oxy_neighbor in carbon_neighbor.GetNeighbors():
                                          if oxy_neighbor.GetAtomicNum() == 16 and oxy_neighbor == s_atom:
                                            is_sulfone_or_similar = True
                                            break
                                      if is_sulfone_or_similar:
                                          break
                              if not is_sulfone_or_similar:
                                return True, "-SR replaces OH- group"
    
    return False, "Not a thiosugar"