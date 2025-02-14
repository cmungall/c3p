"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1]") # check for -COOH and -COO-
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the carbon atom in carboxylic acid group
    carboxylic_carbon = mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1]"))[0]
    carboxylic_atom = mol.GetAtomWithIdx(carboxylic_carbon)


    # Check for direct link to C atom
    linked_to_carbon = False
    for neighbor in carboxylic_atom.GetNeighbors():
      if neighbor.GetAtomicNum() == 6:
        linked_to_carbon = True
        start_carbon = neighbor.GetIdx()
        break
    if not linked_to_carbon:
      return False, "Carboxylic acid group not linked to carbon chain"


    # Count carbons in the chain starting from start_carbon
    carbon_count = 0
    current_atom_idx = start_carbon
    visited_atoms = set()
    
    while True:
        atom = mol.GetAtomWithIdx(current_atom_idx)
        if atom.GetAtomicNum() != 6:
            break  # Stop if not carbon
        carbon_count += 1
        visited_atoms.add(current_atom_idx)
        
        next_carbon = None
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms:
                next_carbon = neighbor.GetIdx()
                break
        if next_carbon is None:
            break
        current_atom_idx = next_carbon
    

    # Check if the chain length is within the range of 13 to 22 carbons.
    if carbon_count < 13 or carbon_count > 22:
       return False, f"Chain length is {carbon_count}, not between 13 and 22 carbons"


    # Check for any non-C/O atoms in the chain (not necessary but can be added)
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8 and atom.GetIdx() not in visited_atoms and atom.GetAtomicNum() != 1:
        return False, "Non-carbon and non-oxygen atoms outside of acid and not H found"

    # Count rotatable bonds - can be used as a sanity check
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
      return False, "Chain not long enough to be a fatty acid"

    return True, f"Molecule is a long-chain fatty acid with a chain length of {carbon_count} carbons."