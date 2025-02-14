"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:35617 tertiary amine
A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of nitrogen atoms
    n_atoms = mol.GetNumAtoms()
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count == 0:
        return False, "No nitrogen atoms present"

    # Check for tertiary nitrogen(s)
    is_tertiary = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            neighbors = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if sum(neighbor != 1 for neighbor in neighbors) == 3:  # 3 non-hydrogen neighbors
                is_tertiary = True
                break

    if not is_tertiary:
        return False, "No tertiary nitrogen atom found"

    # Check for hydrocarbyl groups (alkyl or aryl)
    has_hydrocarbyl = False
    for bond in mol.GetBonds():
        atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if atom1.GetAtomicNum() == 7 or atom2.GetAtomicNum() == 7:  # Nitrogen atom
            if any(neighbor.GetIsAromatic() for neighbor in [atom1, atom2]):
                has_hydrocarbyl = True
                break
            neighbors = [mol.GetAtomWithIdx(neighbor).GetAtomicNum() for neighbor in atom1.GetNeighbors() + atom2.GetNeighbors()]
            if any(neighbor == 6 for neighbor in neighbors):  # Carbon neighbor
                has_hydrocarbyl = True
                break

    if not has_hydrocarbyl:
        return False, "No hydrocarbyl groups attached to nitrogen"

    return True, "Contains a tertiary nitrogen atom with hydrocarbyl groups attached"