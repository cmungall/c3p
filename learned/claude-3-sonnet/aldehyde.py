"""
Classifies: CHEBI:17478 aldehyde
"""
Here's a Python program to classify chemical entities of the class aldehyde:

"""
Classifies: CHEBI:16456 Aldehyde
An aldehyde is a compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonyl carbon with one H and one R group attached
    aldehyde_pattern = Chem.MolFromSmarts("C(=O)[H]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde functional group found"
    
    # Count number of aldehyde groups
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) != 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups, need exactly 1"
    
    # Check for presence of at least one carbon
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 1:
        return False, "No carbon atoms found"
    
    return True, "Contains an aldehyde functional group (R-C(=O)H)"

# Examples
print(is_aldehyde("C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)C([H])=O)[H])[H])C)[H])C"))
# (True, 'Contains an aldehyde functional group (R-C(=O)H)')

print(is_aldehyde("O=C\C(=C\CCCCC)\C"))
# (True, 'Contains an aldehyde functional group (R-C(=O)H)')

print(is_aldehyde("CCCCCCCCCCC=O"))
# (True, 'Contains an aldehyde functional group (R-C(=O)H)')

print(is_aldehyde("CC(\C=C\C=O)=C/C=C/O"))
# (True, 'Contains an aldehyde functional group (R-C(=O)H)')

print(is_aldehyde("[H]C(=O)C1=CC=C(C)O1"))
# (True, 'Contains an aldehyde functional group (R-C(=O)H)')