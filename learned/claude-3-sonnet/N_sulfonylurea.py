"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas contain a urea group where one nitrogen is connected to a sulfonyl group.
    General structure: R-SO2-NH-C(=O)-NH-R'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core N-sulfonylurea pattern
    # [#7] represents any nitrogen
    # The pattern captures R-SO2-NH-C(=O)-NH-R' including variations
    nsulf_pattern = Chem.MolFromSmarts('[#7][#6](=[#8])[#7]S(=O)(=O)[#6,#7,#8]')
    
    # Alternative pattern for reversed arrangement
    nsulf_pattern2 = Chem.MolFromSmarts('[#6,#7,#8]S(=O)(=O)[#7][#6](=[#8])[#7]')
    
    # Pattern for cyclic variations
    cyclic_pattern = Chem.MolFromSmarts('[#7]1[#6](=[#8])[#7]S(=O)(=O)[#6,#7,#8]1')

    if not (mol.HasSubstructMatch(nsulf_pattern) or 
            mol.HasSubstructMatch(nsulf_pattern2) or
            mol.HasSubstructMatch(cyclic_pattern)):
        return False, "Missing required N-sulfonylurea substructure"

    # Exclude invalid cases
    
    # Pattern for direct S-N-N connection without carbonyl
    invalid_pattern1 = Chem.MolFromSmarts('[N]S(=O)(=O)[N]')
    # Pattern for sulfamate derivatives
    invalid_pattern2 = Chem.MolFromSmarts('[O]S(=O)(=O)[N]')
    
    if mol.HasSubstructMatch(invalid_pattern1):
        # Verify that if S-N-N exists, it's part of proper N-sulfonylurea
        for match in mol.GetSubstructMatches(invalid_pattern1):
            n1, s, n2 = match[0], match[1], match[2]
            n1_atom = mol.GetAtomWithIdx(n1)
            n2_atom = mol.GetAtomWithIdx(n2)
            
            # Check if either nitrogen is part of a carbonyl group
            has_carbonyl = False
            for n_atom in [n1_atom, n2_atom]:
                for neighbor in n_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        for c_neighbor in neighbor.GetNeighbors():
                            if c_neighbor.GetSymbol() == 'O' and c_neighbor.GetTotalNumHs() == 0:
                                has_carbonyl = True
                                break
            if not has_carbonyl:
                return False, "Invalid S-N-N arrangement without proper carbonyl group"

    # Additional validation of connectivity
    for match in mol.GetSubstructMatches(nsulf_pattern):
        n1, c, o, n2, s = match[:5]
        # Verify oxidation states and connectivity
        if (mol.GetAtomWithIdx(o).GetTotalNumHs() == 0 and  # Carbonyl oxygen
            mol.GetAtomWithIdx(s).GetTotalValence() == 6):   # Sulfur(VI)
            return True, "Contains N-sulfonylurea group with validated connectivity"
            
    for match in mol.GetSubstructMatches(nsulf_pattern2):
        if len(match) >= 5:  # Ensure we have enough atoms in match
            return True, "Contains N-sulfonylurea group with validated connectivity"
            
    for match in mol.GetSubstructMatches(cyclic_pattern):
        if len(match) >= 5:  # Ensure we have enough atoms in match
            return True, "Contains N-sulfonylurea group with validated connectivity"

    return False, "Structure lacks proper N-sulfonylurea connectivity"