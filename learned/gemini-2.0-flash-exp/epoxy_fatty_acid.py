"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a fatty acid containing an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide ring pattern (more precise definition of an epoxide)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"


    #Check number of oxygen atoms. Fatty acids usually have very few.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 5:
        return False, "Too many oxygen atoms for a fatty acid."


    #Check number of rings that contain both the epoxide and carboxylic groups
    for ep_match in epoxide_matches:
       ep_atoms = [mol.GetAtomWithIdx(i) for i in ep_match]
       for acid_match in acid_matches:
        acid_atoms = [mol.GetAtomWithIdx(i) for i in acid_match]

        for ep_atom in ep_atoms:
            for acid_atom in acid_atoms:
                if ep_atom.GetSymbol() == "C" and acid_atom.GetSymbol() == "C":
                    
                    #if there is a ring including the epoxide carbon and the carboxylic carbon, then its not a fatty acid
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                       if ep_atom.GetIdx() in ring and acid_atom.GetIdx() in ring:
                            return False, "Epoxide and carboxylic acid are part of a ring"



    # Check if the epoxide and the acid are linked by a carbon chain
    # We look for a chain of at least 4 carbons connected to the epoxide and the acid group
    # Create SMARTS pattern for a chain linking the epoxide to the acid. We start from the carbons of the epoxide
    for ep_match in epoxide_matches:
      ep_atoms = [mol.GetAtomWithIdx(i) for i in ep_match]
      for acid_match in acid_matches:
        acid_atoms = [mol.GetAtomWithIdx(i) for i in acid_match]
        for ep_atom in ep_atoms:
            for acid_atom in acid_atoms:
                 if ep_atom.GetSymbol() == "C" and acid_atom.GetSymbol() == "C": # both ends must be carbon
                    
                   # Check if there is a path between the epoxide carbon and the acid carbon
                   #Now we will check if there is a carbon chain of at least 3 carbons
                   
                   #search for a general path between them. We use recursive method.
                   def find_path_recursive(current_atom, target_atom, visited, path):
                    
                        visited.add(current_atom.GetIdx())
                        path.append(current_atom)

                        if current_atom.GetIdx() == target_atom.GetIdx():
                             return path
                        
                        for neighbor in current_atom.GetNeighbors():
                            if neighbor.GetIdx() not in visited and neighbor.GetSymbol() == "C":
                                result = find_path_recursive(neighbor, target_atom, visited.copy(), path.copy())
                                if result:
                                    return result

                        return None
                   
                   path = find_path_recursive(ep_atom, acid_atom, set(), [])

                   if path is None:
                        continue #no path between this epoxide/acid pair, try the next


                   if len(path) < 5: # at least a chain of 3 carbons between them. the 2 end atoms are C too
                        continue

                   all_carbons = True
                   for atom_in_path in path[1:-1]:
                       if atom_in_path.GetSymbol() != "C":
                            all_carbons = False
                            break
                   if not all_carbons:
                        continue

                   
                   # Count rotatable bonds to verify long chains, more flexible way of counting carbons
                   n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
                   if n_rotatable < 6:
                        return False, "Chain too short for a fatty acid"

                   #Check the number of carbons, must be greater than 12.
                   c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                   if c_count < 12:
                        return False, "Too few carbons for a fatty acid"
                    
                   # check the number of double bonds
                   double_bond_count = 0
                   for bond in mol.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            double_bond_count += 1
                   if double_bond_count > 6:
                        return False, "Too many double bonds for a fatty acid"


                   return True, "Contains an epoxide ring and a fatty acid chain"
    
    # if we reach here, it means no path was found
    return False, "No fatty acid chain linking the epoxide and the carboxylic acid group"