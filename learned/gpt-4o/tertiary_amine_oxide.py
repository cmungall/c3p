"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is characterized by a positively charged nitrogen atom 
    bonded to three organic groups, with one of its lone pairs forming a bond with an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nitrogens with three organic groups and one oxygen with negative charge
    n_atom_id = -1
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 4:  # Expecting three organic groups + oxygen with bond
                oxygen_count = sum(1 for n in neighbors if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1)
                if oxygen_count == 1:
                    c_count = sum(1 for n in neighbors if n.GetAtomicNum() in [6] and n.GetFormalCharge() == 0)
                    if c_count == 3:  # Three organic groups
                        n_atom_id = atom.GetIdx()
                        break

    if n_atom_id != -1:
        return True, "Contains tertiary amine oxide structure with correct arrangement"
        
    return False, "Does not contain the correct tertiary amine oxide structure"

# Example usage
test_smiles = "C[N+](C)([O-])C"  # Trimethylamine N-oxide
result, reason = is_tertiary_amine_oxide(test_smiles)
print(f"Result: {result}, Reason: {reason}")