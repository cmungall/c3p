"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine consists of an L-carnitine core with an acyl group attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the L-carnitine core with correct stereochemistry using SMARTS
    # The key here is the @ which designates the chirality
    # The ester oxygen (O) and quaternary nitrogen (N) are linked to the chiral carbon (*)
    # The 3 other carbons (C) have at least one hydrogen.
    carnitine_core_smarts = "[C@H]([OX2])(CC(=O)[O-])C[N+](C)(C)C"
    carnitine_core_pattern = Chem.MolFromSmarts(carnitine_core_smarts)

    if not mol.HasSubstructMatch(carnitine_core_pattern):
        return False, "L-carnitine core not found with the correct stereochemistry"
    
    # Check for an ester group attached to the oxygen of the core
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester bond found"
    
    # check that the ester oxygen atom is indeed the oxygen at the carnitine chiral carbon
    matches = mol.GetSubstructMatches(carnitine_core_pattern)
    found_correct_oxygen = False
    for match in matches:
        chiral_carbon_idx = match[0]
        chiral_carbon = mol.GetAtomWithIdx(chiral_carbon_idx)
        oxygen_idx = -1 # will store the index of the oxygen linked to the chiral carbon
        for neighbor in chiral_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                oxygen_idx = neighbor.GetIdx()
                break
        if oxygen_idx == -1:
           continue #should not happen in a correctly matched structure
        for ester_match in ester_matches:
           if oxygen_idx == ester_match[0]:
                found_correct_oxygen = True
                break
        if found_correct_oxygen:
           break
    
    if not found_correct_oxygen:
        return False, "Ester bond not attached to the carnitine oxygen"


    # check the number of nitrogens (sanity check, should be 1 quaternary)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1)
    if n_count != 1:
      return False, f"Number of quaternary nitrogens incorrect, should be 1, found {n_count}"
    
    
    return True, "Molecule is an O-acyl-L-carnitine"