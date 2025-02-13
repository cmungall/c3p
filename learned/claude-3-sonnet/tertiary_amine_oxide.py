"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine N-oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine N-oxide based on its SMILES string.
    A tertiary amine N-oxide has a nitrogen atom bonded to three organic groups and 
    an oxygen atom with a formal charge separation (N+-O-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine N-oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N+-O- pattern
    n_oxide_pattern = Chem.MolFromSmarts("[NX4+]-[O-]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "No N-oxide group found"
    
    # Get matches for N-oxide groups
    n_oxide_matches = mol.GetSubstructMatches(n_oxide_pattern)
    
    # Check each N-oxide nitrogen
    for match in n_oxide_matches:
        n_idx = match[0]  # Index of nitrogen atom
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check formal charge of N is +1
        if n_atom.GetFormalCharge() != 1:
            continue
            
        # Count number of carbon neighbors (organic groups)
        carbon_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() 
                             if neighbor.GetAtomicNum() == 6)
        
        # For a tertiary amine oxide, need 3 carbon neighbors
        # (the fourth neighbor is the oxide oxygen)
        if carbon_neighbors == 3:
            # Check that nitrogen has no hydrogens
            if n_atom.GetTotalNumHs() == 0:
                return True, "Found nitrogen with three organic groups and N-oxide"
    
    return False, "No nitrogen found with exactly three organic groups and N-oxide"

def test_examples():
    """Test function with some example SMILES"""
    examples = [
        "C[N+](C)([O-])C",  # trimethylamine N-oxide
        "CN(C)C",  # trimethylamine (not N-oxide)
        "C[NH+](C)[O-]",  # N-oxide but secondary amine
        "C[N+](CC)(CC)[O-]",  # triethylamine N-oxide
    ]
    
    for smi in examples:
        result, reason = is_tertiary_amine_oxide(smi)
        print(f"\nSMILES: {smi}")
        print(f"Is tertiary amine oxide: {result}")
        print(f"Reason: {reason}")

if __name__ == "__main__":
    test_examples()