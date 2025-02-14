"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain C18 fatty acid with two C=C double bonds and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    #check that the molecule does not contain any rings.
    ssr = Chem.GetSSSR(mol)
    if len(ssr) > 0:
        return False, "Molecule contains rings"


    # Check for exactly one carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, need exactly 1"


    # Check for exactly two C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
         return False, f"Found {len(double_bond_matches)} double bonds, need exactly 2"

    # Check for 18 carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) != 18:
        return False, f"Found {len(carbon_atoms)} carbon atoms, need exactly 18"
    
    # Check for triple bonds
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Molecule contains triple bonds"

    # Check that every carbon in the chain is connected to at most 2 other carbon atoms (and that the chain is linear)
    carbon_count = 0
    bond_count = 0
    for atom in carbon_atoms:
        carbon_neighbours = [neighbour for neighbour in atom.GetNeighbors() if neighbour.GetAtomicNum() == 6]
        if len(carbon_neighbours) > 2: # Not a linear chain
            
          if len(carbon_neighbours) == 3:
              # check if central carbon of an allene moiety
              if atom.GetTotalValence() != 4:
                return False, "Not a straight carbon chain"

              double_bonds_to_carbon = [bond for bond in atom.GetBonds() if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE) and (bond.GetBeginAtom() in carbon_atoms or bond.GetEndAtom() in carbon_atoms)]
              if len(double_bonds_to_carbon) != 2:
                return False, "Not a straight carbon chain"
              
              
          else:
             return False, "Not a straight carbon chain"
        carbon_count += 1
        for neighbour in carbon_neighbours:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(),neighbour.GetIdx())
            if bond is not None:
                bond_count += 1
            
    # Every bond was counted twice, because it exists between two atoms
    if bond_count != 2*(carbon_count-1):
        return False, "Not a straight carbon chain"
        
    
    
    # Check that one of the carboxylic acid groups is bonded to a terminal carbon (carbon with only one carbon neighbour)
    for match in acid_matches:
      acid_carbon_index = match[0] # index of carbon in COOH
      acid_carbon = mol.GetAtomWithIdx(acid_carbon_index)
      carbon_neighbours = [neighbour for neighbour in acid_carbon.GetNeighbors() if neighbour.GetAtomicNum() == 6]
      if len(carbon_neighbours) != 1: #check that it has only one carbon neighbour, i.e. it is terminal
           return False, "Carboxylic acid group not attached to a terminal carbon."
      
    return True, "Contains 18 carbons in a straight chain, a carboxylic acid group and two double bonds"