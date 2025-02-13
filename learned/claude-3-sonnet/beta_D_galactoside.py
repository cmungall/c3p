"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
A D-galactoside having beta-configuration at its anomeric centre
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_beta_D_galactoside, reason)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Pattern for beta-D-galactopyranose core
    # [OH1,O][C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O
    # Note: The anomeric oxygen can be part of a glycosidic bond
    galactose_pattern = Chem.MolFromSmarts('[OX2][C@H]1O[C@H](C[OX2])[C@H]([OX2])[C@H]([OX2])[C@H]1[OX2]')
    
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No beta-D-galactose moiety found"
    
    # Get matches for the galactose pattern
    matches = mol.GetSubstructMatches(galactose_pattern)
    
    # Check each match
    for match in matches:
        # Get the anomeric carbon (C1) and its attached oxygen atoms
        c1_idx = match[2]  # Index of C1 in the pattern
        c1 = mol.GetAtomWithIdx(c1_idx)
        
        # Get neighboring oxygen atoms of C1
        o_neighbors = [n for n in c1.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(o_neighbors) != 2:
            continue
            
        # One oxygen should be ring oxygen, other should be glycosidic/substituent
        ring_o = None
        subst_o = None
        for o in o_neighbors:
            if len([n for n in o.GetNeighbors() if n.GetAtomicNum() == 6]) == 2:
                ring_o = o
            else:
                subst_o = o
                
        if ring_o is None or subst_o is None:
            continue
            
        # Verify beta configuration by checking if substituent is below the ring
        # In the given pattern, beta configuration is encoded in the stereochemistry
        return True, "Contains beta-D-galactoside moiety"
        
    return False, "No beta-D-galactoside configuration found"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",  # methyl beta-D-galactoside
        "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)C1=O",  # Not a galactoside
        "CO[C@H]1O[C@H](CO)[C@@H](O)[C@@H](O)[C@H]1O"   # methyl alpha-D-glucoside
    ]
    
    for smi in examples:
        result, reason = is_beta_D_galactoside(smi)
        print(f"SMILES: {smi}")
        print(f"Is beta-D-galactoside: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()