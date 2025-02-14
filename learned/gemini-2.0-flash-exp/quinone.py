"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a compound having a fully conjugated cyclic dione structure, 
    derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring(s)
    if mol.GetRingInfo().NumRings() == 0:
        return False, "Molecule does not contain a ring"

    # Check for two carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 2:
        return False, f"Molecule has {len(carbonyl_matches)} carbonyl groups, needs at least 2"
    
    # Conjugation check and aromatic parent check
    ring_systems = []
    for match in carbonyl_matches:
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            for ring_index in mol.GetAtomRingInfo().AtomRings():
                 if atom_index in ring_index:
                    ring = mol.GetRingInfo().GetAtomIndices(ring_index)
                    if ring not in ring_systems:
                       ring_systems.append(ring)
    
    
    conjugated_ring_found = False
    for ring in ring_systems:
        #get the atoms from the ring
        ring_mol = Chem.Mol(mol)
        remove_atoms_ids = []
        for atom in mol.GetAtoms():
           if not atom.GetIdx() in ring:
             remove_atoms_ids.append(atom.GetIdx())

        ring_mol = Chem.RWMol(mol)
        remove_atoms_ids.sort(reverse=True)
        for atom_id in remove_atoms_ids:
            ring_mol.RemoveAtom(atom_id)

        ring_mol = ring_mol.GetMol()
        
        
        #check if the ring is aromatic
        if ring_mol.GetRingInfo().IsRingAromatic(ring):
              conjugated_ring_found = True
              break
        else:
            #check if the ring has at least 3 double bonds
            double_bonds_count = 0
            for bond in ring_mol.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_bonds_count += 1
            if double_bonds_count >= 3:
                conjugated_ring_found = True
                break
        
    if not conjugated_ring_found:
            return False, "No conjugated ring with carbonyl groups found"
    
    
    # Check parent aromatic
    parent_aromatic = False
    for ring in ring_systems:
       ring_mol = Chem.Mol(mol)
       remove_atoms_ids = []
       for atom in mol.GetAtoms():
         if not atom.GetIdx() in ring:
           remove_atoms_ids.append(atom.GetIdx())
       
       ring_mol = Chem.RWMol(mol)
       remove_atoms_ids.sort(reverse=True)
       for atom_id in remove_atoms_ids:
         ring_mol.RemoveAtom(atom_id)
       ring_mol = ring_mol.GetMol()
        
       
       carbonyls_in_ring = False
       for match in carbonyl_matches:
           for atom_index in match:
             if atom_index in ring:
               carbonyls_in_ring = True
               break
           if carbonyls_in_ring:
             break
       if carbonyls_in_ring:
         
          no_carbonyl_ring_mol = Chem.Mol(ring_mol)
          
          carbonyl_atoms_ids_to_remove = []
          for match in carbonyl_matches:
              for atom_index in match:
                if atom_index in ring:
                    carbonyl_atoms_ids_to_remove.append(atom_index)
          
          no_carbonyl_ring_mol = Chem.RWMol(ring_mol)
          carbonyl_atoms_ids_to_remove.sort(reverse=True)
          for atom_id in carbonyl_atoms_ids_to_remove:
              no_carbonyl_ring_mol.RemoveAtom(atom_id)
          no_carbonyl_ring_mol = no_carbonyl_ring_mol.GetMol()
          
          if len(no_carbonyl_ring_mol.GetAtoms()) > 0:
              
             if no_carbonyl_ring_mol.GetRingInfo().IsRingAromatic(ring):
                 parent_aromatic = True
                 break;
       else:
           if ring_mol.GetRingInfo().IsRingAromatic(ring):
               parent_aromatic = True
               break
           
    if not parent_aromatic:
           return False, "Not derived from an aromatic compound"


    return True, "Molecule has a conjugated cyclic dione structure derived from an aromatic compound"