"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is any terpenoid derived from a sesquiterpene, often modified by rearrangement 
    or methyl group removal, having a C15 carbon skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Contains {c_count} carbons, less than required for sesquiterpenoids"

    if c_count > 15:
        return False, f"Contains {c_count} carbons, more complex than a typical sesquiterpenoid"

    # Count the number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "Lacks ring structures typically present in sesquiterpenoids"

    # Look for common functional groups indicating derivations (e.g., -OH, =O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "Does not contain oxygen, unlikely to be a sesquiterpenoid"

    return True, "Contains characteristics typical of a sesquiterpenoid"

# Example
print(is_sesquiterpenoid('O=C1C=C2C=CC(=O)[C@@]([C@]2(C)C[C@]1(O)C(=C)COC(=O)C)(O)C'))  # Example test