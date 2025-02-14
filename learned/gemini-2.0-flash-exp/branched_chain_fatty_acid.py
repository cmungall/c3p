"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a carboxylic acid group, a long carbon chain and one or more alkyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for branching point and long chain. 
    branch_found = False
    chain_length_ok = False
    
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6: # Check only carbon atoms
            neighbors = [neighbor for neighbor in atom.GetNeighbors()]
            carbon_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 6]
            non_h_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() != 1]
            
            # A branch should be connected to at least 2 carbon neighbors and another non-H
            if len(carbon_neighbors) >= 2 and len(non_h_neighbors) > 2:
                branch_found = True
                
    if not branch_found:
      return False, "No appropriate branching point found"
    
    #check rotatable bonds
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 2:
        return False, "Chain is too short"
    
    # Check molecular weight - fatty acids usually > 90 Da (for molecules like 2-methylbutyric acid)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 90:
         return False, "Molecular weight too low for fatty acid"
    
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 3 or o_count != 2:
        return False, "Molecule is too short"

    return True, "Contains a carboxylic acid group, a long chain and one or more branching points"