"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    A cephalosporin has a beta-lactam ring fused to a 6-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core cephalosporin substructure with a 6 membered ring fused to a beta-lactam ring (and 1 S atom)
    # The ? makes the R group variable, 
    core_pattern = Chem.MolFromSmarts('S1C[C@H]2N(C1)C(=O)[C@@H]2[?][?]')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 6-membered dihydrothiazine and beta-lactam ring not found"


    #Check for carboxylic acid
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(carboxylic_matches) == 0:
        return False, "No carboxylic acid group found"

    #Check for acylamino substituent at position 7
    acylamino = Chem.MolFromSmarts('N[C@H](C=O)[C@H]1[N]2C(=O)C([C@H]2[?]S1)[?][?][?]')
    if not mol.HasSubstructMatch(acylamino):
        return False, "No acylamino group on position 7"

    
    return True, "Contains core cephalosporin structure with a carboxylic acid and acylamino group"