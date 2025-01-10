"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin typically contains a 20-carbon skeleton that includes
    a cyclopentane ring and functional groups such as hydroxyl, ketone, and 
    carboxylic acid attached to the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule resembles a prostaglandin structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Attempt to find general prostaglandin likeness via SMARTS
    # Prostaglandin pattern: Cyclopentane with chains and optional functional groups
    prostaglandin_pattern = Chem.MolFromSmarts('[C@H]1[C@H](C=C1)CC(=O)|CC(=O)|O|CC=C|C=C')
    
    if not mol.HasSubstructMatch(prostaglandin_pattern):
        return False, "Structure does not match typical prostaglandin cyclopentane pattern"

    # Count key oxygens for groups - prostaglandins generally have multiple hydrophilic groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:  # Heuristic: Prostaglandins tend to have at least 3 oxygens (OH, COOH, ketones)
        return False, "Too few oxygen atoms for typical prostaglandin"

    # Check overall carbon count, should typically be 20
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Number of carbons is {c_count}, expected about 20 for prostaglandins"

    # Check for significant stereochemistry - many prostaglandins are chiral
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, "Chirality missing; prostaglandins are chiral molecules"

    return True, "Structure contains features typical of prostaglandins: cyclopentane ring, multiple oxygens, specific stereochemistry"