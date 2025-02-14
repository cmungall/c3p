"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a nucleoside 5'-phosphate, False otherwise, and reason
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    phosphate_pattern = Chem.MolFromSmarts("[CX4]-[OP](=[OX1])([OX1])[OX1]") # Phosphate attached to a carbon
    ribose_pattern = Chem.MolFromSmarts("C1[CH][CH](O)[CH](O)[O]1")  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts("C1[CH][CH](O)[CH](C)[O]1") #Deoxyribose
    purine_pattern = Chem.MolFromSmarts("c1nc2c(n1)ncn2") # Purine
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncc(=O)[nH]1") # Pyrimidine

    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
          return False, "No phosphate group found"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    # Check that at least one of the phosphates is attached to the C5' of the sugar
    phosphate_c_atoms = []
    for match in phosphate_matches:
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "C":
                phosphate_c_atoms.append(atom_index) #store all carbon atoms attached to phosphate


    has_5prime_phosphate = False
    for c_atom_index in phosphate_c_atoms:
        c_atom = mol.GetAtomWithIdx(c_atom_index) # get atom
        neighbors = c_atom.GetNeighbors()
        for neighbor in neighbors:
          if neighbor.GetSymbol() == "O":
            for n_n in neighbor.GetNeighbors():
                if n_n.GetSymbol() == "C":
                    if (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
                      if (len(n_n.GetNeighbors()) == 4): #this is C5
                        has_5prime_phosphate = True
                        break
                if has_5prime_phosphate:
                  break;
          if has_5prime_phosphate:
            break;
        if has_5prime_phosphate:
          break;
    if not has_5prime_phosphate:
      return False, "Phosphate not at 5' position"
    
    # Check for ribose or deoxyribose ring
    if not (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
        return False, "No ribose/deoxyribose ring found"

    # Check for purine or pyrimidine base attached to the sugar
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No purine/pyrimidine base found"
    
    # If all conditions are met, return True
    return True, "Molecule is a nucleoside 5'-phosphate"