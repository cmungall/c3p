"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin must have a penam core structure with additional specific substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define penam core pattern
    # Penam core: a four-membered beta-lactam ring fused to a five-membered thiazolidine ring
    penam_core_pattern = Chem.MolFromSmarts('C1C([N+]2C(=C(C2=O)S1)C)(C)C')
    if not mol.HasSubstructMatch(penam_core_pattern):
        return False, "No penam core structure found"

    # Check for two methyl groups at position 3
    methyl_pattern = Chem.MolFromSmarts('C1(C)C(S1)(C)')
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Two methyl groups at position 3 not found"

    # Check for carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts('C(=O)[O-]')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Carboxylate group at position 3 not found"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Carboxamido group at position 6 not found"
    
    # Ensure the molecule also matches specific stereochemistry
    stereo = rdMolDescriptors.CalcMolStereo(mol)
    if stereo < 2:  # Just a rough check on stereochemistry for validity; refine as needed
        return False, "Stereochemistry does not match penicillin"

    return True, "Molecule matches the penicillin structure"

# Test the function 
smiles_example = '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)C(O)=O'
print(is_penicillin(smiles_example))