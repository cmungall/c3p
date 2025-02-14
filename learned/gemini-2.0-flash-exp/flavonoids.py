"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid contains the C6-C3-C6 core structure, typically as a benzopyran derivative,
    and includes various subgroups such as flavones, flavanones, isoflavones, chalcones, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core C6-C3-C6 (two phenyl rings connected by a 3 carbon chain)
    # using smarts to account for any of the 2 rings being heterocyclic
    core_pattern_1 = Chem.MolFromSmarts("[c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1[C][C][C][c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1")
    core_pattern_2 = Chem.MolFromSmarts("[c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1[C](=[O])[C][C][c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1")
    if not (mol.HasSubstructMatch(core_pattern_1) or mol.HasSubstructMatch(core_pattern_2)):
         # Attempting a more generic C6-C3-C6
        core_pattern_generic = Chem.MolFromSmarts("[c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1[C,c][C,c][C,c][c]1[c,o,n][c,o,n][c,o,n][c,o,n][c,o,n]1")
        if not mol.HasSubstructMatch(core_pattern_generic):
            return False, "No C6-C3-C6 core structure found"


    # Heterocyclic Ring: Look for a common benzopyran-like structure
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2c(cc1)occc2") #simple benzopyran
    chromene_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=[OX1])cc2") # chromene
    if not (mol.HasSubstructMatch(benzopyran_pattern) or mol.HasSubstructMatch(chromene_pattern)):
        return False, "No common benzopyran or chromene-like ring system found"
    
    #Check if there are at least 2 phenyl rings
    phenyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("c1ccccc1")))
    if phenyl_count < 2:
         return False, "Molecule does not have at least 2 phenyl rings."

    # Check for common oxygen substitutions
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
       return False, "Not enough oxygen atoms for flavonoid."
    
    return True, "Likely a flavonoid based on core structure and heterocyclic ring"