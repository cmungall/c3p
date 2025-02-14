"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    A omega-hydroxy fatty acid is defined as a straight-chain fatty acid with a carboxyl group at position 1
    and a hydroxyl group at the omega (last) position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (C(=O)O or C(O)=O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
         carboxyl_pattern = Chem.MolFromSmarts("[C;X3]([OH1])[OH0]")  # Handles tautomer
         carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
         if not carboxyl_matches:
            return False, "No carboxyl group found"

    # Get the index of the carboxyl carbon
    carbonyl_carbon_idx = carboxyl_matches[0][0]


    # Get atoms
    atoms = mol.GetAtoms()
    
    # Check if the molecule has at least 3 carbons. If not, it cannot be a fatty acid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
       return False, "Too few carbon atoms"

    # Start from the carboxyl carbon and traverse the chain, getting the last carbon of the chain
    current_atom_idx = carbonyl_carbon_idx
    previous_atom_idx = -1 #Initialize to not be any real atom index

    # Traverse until a terminal carbon is reached
    while True:
        current_atom = atoms[current_atom_idx]
        neighbors = [n.GetIdx() for n in current_atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != previous_atom_idx] #Get only carbon neighbors
        
        if len(neighbors) == 0:
          omega_carbon_idx = current_atom_idx
          break  # Terminal carbon reached
        elif len(neighbors) == 1:
           previous_atom_idx = current_atom_idx
           current_atom_idx = neighbors[0]
        else: #Branching, therefore not a straight chain fatty acid
            return False, "Not a straight-chain fatty acid"

    #Check if an -OH group is present at the omega position
    omega_carbon = atoms[omega_carbon_idx]
    
    for neighbor in omega_carbon.GetNeighbors():
       if neighbor.GetAtomicNum() == 8 and (neighbor.GetTotalNumHs() == 1):
           return True, "Omega hydroxyl group found"

    return False, "No omega-hydroxyl group found"