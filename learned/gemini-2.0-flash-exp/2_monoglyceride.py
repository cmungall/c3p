"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol with one fatty acid attached to the central (2nd) carbon via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])[CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
       return False, "No glycerol backbone found"
    
    # Check for *one* ester group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
       return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check if hydroxyl groups are free (not esterified) in 1 and 3 position.
    free_hydroxyl_pattern = Chem.MolFromSmarts("[CH2X4]([OX2])")
    hydroxyl_matches = mol.GetSubstructMatches(free_hydroxyl_pattern)

    # Get the atoms of the glycerol backbone
    glycerol_atoms = list(glycerol_matches[0])

    # Check that the ester is on the middle C
    ester_connected = False
    for match in ester_matches:
        ester_o_atom_index = match[0]
        ester_o_atom = mol.GetAtomWithIdx(ester_o_atom_index)
        
        for neighbor in ester_o_atom.GetNeighbors():
            if neighbor.GetIdx() == glycerol_atoms[1]:  # glycerol atoms: 0-1-2, so 1 is the middle carbon
                ester_connected = True
                break
        if ester_connected:
            break

    if not ester_connected:
        return False, "Ester group not connected to the central carbon of the glycerol"


    # Verify that the other 2 hydroxyls are free
    free_OH_count = 0
    for match in hydroxyl_matches:
      if match[0] == glycerol_atoms[0] or match[0] == glycerol_atoms[2]:
        free_OH_count +=1
    
    if free_OH_count != 2:
        return False, "Did not find two free hydroxyl groups at the 1 and 3 positions of the glycerol"


    # Check for fatty acid chain connected to the ester
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # at least 4 C in a row
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be a fatty acid."
    
    # Check if all oxygen are part of OH or ester linkages
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 4:
      return False, "Must have exactly four oxygens (2 hydroxyl, 1 ester)"


    return True, "Contains a glycerol backbone with one fatty acid attached at the 2nd carbon via ester bond"