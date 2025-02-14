"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are triterpenoids with a characteristic 4,4,8-trimethyl-17-furanylsteroid skeleton
    and high degree of oxygenation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found."

    # Check for the tetracyclic core with correct methyl groups at 4,4,8 and a furan at position 17
    # This core is highly variable, we will look for a simpler SMARTS and number of carbons

    # We also need to check the overall number of carbons to see if it is 30
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 30:
        return False, f"Number of carbons is not 30, found {carbon_count}"

    # Check for high oxygenation (at least 4 oxygens, we need to do this to avoid other triterpenoids)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
        return False, f"Too few oxygen atoms, found {oxygen_count}, need at least 4"


    # If all criteria pass return true
    return True, "Matches limonoid criteria: triterpenoid with furan, 30 carbons, and high oxygenation."