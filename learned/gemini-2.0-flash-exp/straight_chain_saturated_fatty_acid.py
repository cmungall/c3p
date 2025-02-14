"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Reason for the classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Verify Carboxylic Acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for saturation
    if not all(bond.GetBondType() == Chem.BondType.SINGLE for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6):
        return False, "Not a saturated molecule"
    
    # Check for allowed elements (C, H, O and D)
    allowed_elements = [1, 6, 8, 2]  # H, C, O, D
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Molecule contains disallowed elements"

    # Check for branching, by looking for carbons with more than two carbon neighbors.
    # The carbonyl carbon of the carboxylic group is excluded from the branching check
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetSymbol() == 'C': #carbon
            if not atom.HasProp('carbonyl') or not atom.GetProp('carbonyl') == 'True':
                carbon_neighbors = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                       carbon_neighbors += 1
                if carbon_neighbors > 2 :
                    return False, "Branched chain detected"
                
    # Check number of carbon atoms (must be at least 4)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbon atoms"

    # Check for hydroxy groups, allowing 0 or 1 of them attached to a carbon in the chain.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    n_hydroxy = len(hydroxy_matches)
    if n_hydroxy > 1:
        return False, "Too many hydroxy groups"
    
    if n_hydroxy == 1:
       for match in hydroxy_matches:
           hydroxy_atom = mol.GetAtomWithIdx(match[0])
           
           #check is oxygen is connected to a carbon
           oxygen_neighbors = hydroxy_atom.GetNeighbors()
           if not any(neighbor.GetAtomicNum() == 6 for neighbor in oxygen_neighbors):
               return False, "Hydroxy group not connected to carbon chain"

           #check that the carbon attached to the hydroxy group is not the carbonyl group
           for neighbor in oxygen_neighbors:
               if neighbor.GetAtomicNum() == 6 and neighbor.HasProp('carbonyl') and neighbor.GetProp('carbonyl') == 'True':
                    return False, "Hydroxy group attached to carbonyl carbon"
    
    # Check for other functional groups.
    # First, label all carbonyl carbons in carboxyl groups
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    for match in carbonyl_matches:
        mol.GetAtomWithIdx(match[0]).SetProp('carbonyl', 'True')


    # Then, check for any other non-hydrogen atoms which are not carbons, oxygens in the carboxyl group or in the optional hydroxy group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 6: # not hydrogen or carbon
           if atom.GetAtomicNum() == 8: # if oxygen
              if not any(neighbor.GetAtomicNum() == 6 and neighbor.HasProp('carbonyl') and neighbor.GetProp('carbonyl') == 'True' for neighbor in atom.GetNeighbors()): #is not the oxygen in a carboxyl group
                 if not any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()): #is not the oxygen in a hydroxy group
                    return False, "Other functional groups detected"
           else:
               return False, "Other functional groups detected"

    return True, "Straight-chain saturated fatty acid"