"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    A sphingomyelin d18:1 is any sphingomyelin having sphingosine as the sphingoid component.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for phosphocholine group
    phosphocholine_smarts = '[O]-[P](=O)([O-])-[O][CH2][CH2][N+](C)(C)C'
    phosphocholine_pattern = Chem.MolFromSmarts(phosphocholine_smarts)
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
    
    # Check for amide bond connected to fatty acyl chain
    amide_smarts = '[NX3][C](=O)[C]'
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond to fatty acyl chain found"
    
    # Identify sphingosine backbone
    # Find the amide nitrogen atom
    amide_nitrogen_idx = None
    for match in amide_matches:
        amide_nitrogen_idx = match[0]
        break  # Assume first match is the sphingoid amide
    if amide_nitrogen_idx is None:
        return False, "Amide nitrogen not found"
    
    # Find the oxygen atom connected to phosphocholine group
    phosphocholine_matches = mol.GetSubstructMatches(phosphocholine_pattern)
    phospho_oxygen_idx = None
    for match in phosphocholine_matches:
        phospho_oxygen_idx = match[0]  # The oxygen atom connected to phosphate
        break
    if phospho_oxygen_idx is None:
        return False, "Phosphocholine attachment oxygen not found"
    
    # Find the path between amide nitrogen and phosphocholine oxygen
    path = Chem.rdmolops.GetShortestPath(mol, amide_nitrogen_idx, phospho_oxygen_idx)
    # Remove the amide bond (amide nitrogen and its attached carbonyl carbon)
    path = list(path)
    if len(path) == 0:
        return False, "No path between amide nitrogen and phosphocholine oxygen"
    
    # Analyze the backbone chain
    backbone = mol.GetSubstructMatch(Chem.MolFragmentToSmarts(mol, atomsToUse=path))
    if not backbone:
        return False, "Sphingosine backbone not found"
    
    # Count the number of carbons in the backbone
    carbon_count = 0
    double_bond_count = 0
    hydroxyl_count = 0
    amino_count = 0
    for idx in path:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
        elif atom.GetAtomicNum() == 8:
            # Check if it's a hydroxyl group attached to carbon
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if 6 in neighbors:
                hydroxyl_count += 1
        elif atom.GetAtomicNum() == 7:
            amino_count += 1
    
    # Check for double bonds in the backbone
    for i in range(len(path)-1):
        bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
    
    if carbon_count != 18:
        return False, f"Backbone has {carbon_count} carbons, expected 18"
    if double_bond_count != 1:
        return False, f"Backbone has {double_bond_count} double bonds, expected 1"
    if hydroxyl_count < 2:
        return False, f"Backbone has {hydroxyl_count} hydroxyl groups, expected at least 2"
    if amino_count < 1:
        return False, "Backbone missing amino group"
    
    return True, "Molecule is a sphingomyelin d18:1 with sphingosine backbone"

__metadata__ = {
    'chemical_class': {
        'name': 'sphingomyelin d18:1',
        'definition': 'Any sphingomyelin having sphingosine as the sphingoid component.'
    },
    'message': None,
    'success': True
}