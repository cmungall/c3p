"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid with one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for at least one hydroxyl group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(alcohol_pattern):
         return False, "No hydroxy group found"

    # Get the carbon of the carboxylic acid
    carboxyl_carbon_matches = mol.GetSubstructMatches(acid_pattern)
    carboxyl_carbon = carboxyl_carbon_matches[0][0]

    # Check for a long aliphatic carbon chain connected to the carboxylic acid
    def check_chain_length(start_atom):
      """Checks the length of a carbon chain starting from a given atom."""
      chain_length = 0
      current_atom = start_atom
      visited_atoms = {current_atom}
      
      while True:
        neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(current_atom).GetNeighbors() if n.GetAtomicNum() == 6]
        
        next_atom = None
        for neighbor in neighbors:
            if neighbor not in visited_atoms:
               next_atom = neighbor
               break

        if next_atom is None:
           break
        
        if len([n.GetIdx() for n in mol.GetAtomWithIdx(next_atom).GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
           return False # Found branching, so this isn't the main chain.

        chain_length += 1
        visited_atoms.add(next_atom)
        current_atom = next_atom
          
      return chain_length >= 6 # Chain must be at least 6 carbons.
    
    # Check for chain in both directions.
    chain_carbons = [n.GetIdx() for n in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors() if n.GetAtomicNum() == 6]

    long_chain_found = False
    for carbon in chain_carbons:
        if check_chain_length(carbon):
            long_chain_found = True
            break

    if not long_chain_found:
        return False, "No long carbon chain attached to carboxylic group found"

    # Check for at least ONE hydroxy group on the fatty acid chain
    hydroxy_on_chain = False
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 6:
        for neighbor in atom.GetNeighbors():
          if neighbor.GetSymbol() == 'O' and neighbor.GetTotalValence() == 2:
            # Check if the hydroxyl-bearing carbon is part of the long chain
            
            chain_atom = atom.GetIdx()
            
            def is_on_chain(start_atom, hydroxy_bearing_carbon):
              """Determines if a given carbon is part of the chain, based on the starting atom"""
              current_atom = start_atom
              visited_atoms = {current_atom}
              
              while True:
                if current_atom == hydroxy_bearing_carbon:
                    return True

                neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(current_atom).GetNeighbors() if n.GetAtomicNum() == 6]
                next_atom = None
                for neighbor in neighbors:
                    if neighbor not in visited_atoms:
                      next_atom = neighbor
                      break
                
                if next_atom is None:
                    return False
                
                if len([n.GetIdx() for n in mol.GetAtomWithIdx(next_atom).GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
                   return False # Found branching, so this isn't the main chain.
                
                visited_atoms.add(next_atom)
                current_atom = next_atom
              
              
            for carbon in chain_carbons:
              if is_on_chain(carbon, chain_atom):
                  hydroxy_on_chain = True
                  break
        if hydroxy_on_chain:
          break

    if not hydroxy_on_chain:
        return False, "Hydroxy is not on the fatty acid chain"
            
    # Check the number of carbons:
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Too few carbons to be a fatty acid"

    return True, "Contains a carboxylic acid group, at least one hydroxy group on the fatty acid chain and a long carbon chain"