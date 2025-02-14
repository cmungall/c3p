"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a 4'-hydroxyflavanone, False otherwise.
         str:  Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flavanone core pattern (including numbering from the examples)
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccc2C(=O)CC(c3ccccc3)Oc2c1")

    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Molecule does not contain a flavanone core"

    # Pattern to find the 4'-hydroxy substitution, including where the R group is not carbon.
    hydroxy_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")

    # Find substructures in the molecule that match phenyl rings.
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
    
    found_4_hydroxy = False
    for match in phenyl_matches:
      # Get the phenyl ring atoms
      phenyl_ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
      
      # Check the degree of substitution of each atom
      degrees = {atom.GetIdx():atom.GetDegree() for atom in phenyl_ring_atoms}
      
      # Look for an atom with degree 3 with a hydroxyl group bonded directly
      for atom_idx, degree in degrees.items():
        if degree == 3:
          atom = mol.GetAtomWithIdx(atom_idx)
          for neighbor in atom.GetNeighbors():
              if neighbor.GetAtomicNum() == 8 and len(neighbor.GetNeighbors()) == 1 and (neighbor.GetNeighbors()[0].GetIdx() == atom_idx): #check neighbour is -O and if O is bonded only to the current atom

                
                # Check which position in the ring this atom is. This is done using the SMARTS pattern with dummy atoms
                # To do that, we generate a substructure to match based on the atom position we found.
                
                smarts = "[cH1]1~[cH1]~[cH1]~[cH1]([OH1])~[cH1]~[cH1]1"
                dummy_pattern = Chem.MolFromSmarts(smarts)


                # now generate a list of all atoms we're interested in from the match on the phenyl ring.
                match_atoms = [phenyl_ring_atoms[0].GetIdx(), phenyl_ring_atoms[1].GetIdx(), phenyl_ring_atoms[2].GetIdx(), atom_idx, phenyl_ring_atoms[4].GetIdx(), phenyl_ring_atoms[5].GetIdx()]


                # Now, use this to create a new molecule, which we will compare against the generated pattern
                submol = Chem.EditableMol(mol)
                for atom in mol.GetAtoms():
                    if atom.GetIdx() not in match_atoms:
                        submol.RemoveAtom(atom.GetIdx())
                submol = submol.GetMol()
                                        
                if submol.HasSubstructMatch(dummy_pattern):
                    found_4_hydroxy = True
                    break
          if found_4_hydroxy:
              break

    if found_4_hydroxy:
      return True, "Molecule contains a flavanone core with a hydroxy substituent at the 4' position of the phenyl ring"
    else:
        return False, "Molecule does not have a hydroxyl substituent at the 4' position of the phenyl ring"