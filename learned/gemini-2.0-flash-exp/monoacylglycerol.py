"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with one fatty acid chain attached via ester bonds.
    The other two OH groups may or not be linked to other things that are not fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Find glycerol backbone carbons.
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No glycerol backbone found"

    glycerol_carbon_indices = list(glycerol_match)
    glycerol_carbons = [mol.GetAtomWithIdx(i) for i in glycerol_carbon_indices]

    # Check that each glycerol carbon is bound to at least one oxygen and that there are 3 such oxygens.
    glycerol_oxygen_indices = []
    for carbon in glycerol_carbons:
        o_count = 0
        for neighbor in carbon.GetNeighbors():
           if neighbor.GetSymbol() == "O":
              glycerol_oxygen_indices.append(neighbor.GetIdx())
              o_count += 1

    if len(glycerol_oxygen_indices) < 3:
       return False, "Glycerol backbone not properly substituted with 3 oxygens"

    # 2. Check that one and only one O is part of an ester bond (only considering those bonded to glycerol carbons)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_count = 0
    ester_o_index = -1
    
    for o_index in glycerol_oxygen_indices:
        o_atom = mol.GetAtomWithIdx(o_index)
        for match in mol.GetSubstructMatches(ester_pattern):
            if o_index == match[0]:
               ester_count += 1
               ester_o_index = o_index
               break
        if ester_count > 1:
            break

    if ester_count != 1:
         return False, f"Found {ester_count} ester groups attached to glycerol, need exactly 1"


    # 3. Check for the presence of a fatty acid chain attached to the ester oxygen of the glycerol backbone
    if ester_o_index == -1:
        return False, "Could not identify the ester oxygen"
    
    ester_o_atom = mol.GetAtomWithIdx(ester_o_index)
    fatty_acid_found = False
    for neighbor in ester_o_atom.GetNeighbors():
        if neighbor.GetSymbol() == "C" and neighbor.GetHybridization() == Chem.HybridizationType.SP2:  # Carbonyl carbon
            for carbonyl_neighbor in neighbor.GetNeighbors():
              if carbonyl_neighbor.GetSymbol() == "O":
                 for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetSymbol() == "C": #start of a fatty acid chain
                     
                      current_atom = next_neighbor
                      visited_atoms = {current_atom.GetIdx()}
                      chain_length = 0
                      while True:
                          is_chain_atom = True
                          chain_length += 1
                          
                          next_atom = None
                          for next_neighbor in current_atom.GetNeighbors():
                              if next_neighbor.GetSymbol() == "C" and next_neighbor.GetIdx() not in visited_atoms:
                                  next_atom = next_neighbor
                                  break
                              elif next_neighbor.GetSymbol() != 'H' and next_neighbor.GetSymbol() != 'C':
                                  is_chain_atom=False
                                  break
          
                          if not is_chain_atom:
                              break
                          if next_atom is None:
                              break
                          else:
                              visited_atoms.add(next_atom.GetIdx())
                              current_atom = next_atom
                      if chain_length > 3:
                            fatty_acid_found = True
                            break
            if fatty_acid_found:
                break


    if not fatty_acid_found:
        return False, "Missing fatty acid chain"


    # 4. Check number of rotatable bonds (less strict now)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1:
        return False, "Chain too short to be fatty acid"
    

     # 5. Check the number of carbons and oxygens (also less strict)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 4:
        return False, "Too few carbons for monoacylglycerol"
    if o_count < 3:
        return False, "Must have at least 3 oxygens (glycerol backbone and one ester group)"

    return True, "Contains glycerol backbone with 1 fatty acid chain attached via ester bonds"