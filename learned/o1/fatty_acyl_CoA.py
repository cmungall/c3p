"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is an acyl-CoA that results from the formal condensation 
    of the thiol group of coenzyme A with the carboxy group of any fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    # Check for thioester linkage
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Define SMARTS pattern for adenosine moiety in coenzyme A
    adenosine_pattern = Chem.MolFromSmarts("n1c2ncnc2n(c1)[C@H]3O[C@H](CO)[C@@H](O)[C@H]3O")
    if not mol.HasSubstructMatch(adenosine_pattern):
        return False, "No adenosine moiety found, not coenzyme A"
    
    # Define SMARTS pattern for pantetheine unit in coenzyme A
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found, not coenzyme A"

    # Check for fatty acyl chain attached via thioester linkage
    # Look for long carbon chains attached to the carbonyl of thioester
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)SCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No fatty acyl chain attached via thioester linkage found"

    # Optional: Check length of fatty acyl chain (typically >4 carbons)
    # Get the carbon count in the fatty acyl chain
    fatty_acyl_chain = Chem.MolFromSmarts("C(=O)[CX4]")
    matches = mol.GetSubstructMatches(fatty_acyl_chain)
    fatty_chain_lengths = []
    for match in matches:
        carbon = mol.GetAtomWithIdx(match[1])
        chain_length = 1  # Start with 1 carbon
        # Traverse the carbon chain
        atoms_visited = set()
        atoms_to_visit = [carbon.GetIdx()]
        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in atoms_visited:
                continue
            atoms_visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in atoms_visited:
                        atoms_to_visit.append(neighbor.GetIdx())
        fatty_chain_lengths.append(chain_length)
    if not fatty_chain_lengths or max(fatty_chain_lengths) < 5:
        return False, "Fatty acyl chain is too short"

    return True, "Contains coenzyme A moiety with fatty acyl chain attached via thioester linkage"