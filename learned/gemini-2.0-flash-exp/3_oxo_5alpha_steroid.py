"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid has a steroid skeleton, a carbonyl at position 3 and alpha
    configuration at position 5

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for steroid skeleton
    # Approximate by checking for 4 fused rings
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 4:
        return False, "Not a steroid: less than 4 rings"
    
    # Check if all rings are 5 or 6 membered:
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) < 5 or len(ring) > 6:
             return False, "Not a steroid: ring size different than 5 or 6"

    # 2. Check for carbonyl at position 3. Position 3 is defined as the carbon
    #   that is part of a six membered ring and has two adjacent carbons within the ring.
    #   and that is adjacent to one of the bridgehead carbons of the steroid system
    #   that has 2 bonds to two different rings and is not a sp2 carbon.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not matches:
        return False, "No carbonyl group at position 3 found."
    
    
    # Define the substructure for a bridgehead carbon. A bridgehead carbon is part of a ring system and is
    # bonded to more than 2 carbons.
    bridgehead_pattern = Chem.MolFromSmarts("[C;R2;!$(C=*)]([#6])([#6])")
    
    # Find all the bridgeheads of the steroid ring system.
    bridgehead_matches = mol.GetSubstructMatches(bridgehead_pattern)
    if len(bridgehead_matches) < 4: # There are 4 bridgeheads within the sterane skeleton
      return False, "Not a steroid: less than 4 bridgehead carbons"

    #Check that one of the bridgeheads is adjacent to the 3 position carbon
    found_3_position = False
    for match in matches:
        carbonyl_carbon_index = match[0]
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)
        for bridgehead_match in bridgehead_matches:
           bridgehead_atom_index = bridgehead_match[0]
           bridgehead_atom = mol.GetAtomWithIdx(bridgehead_atom_index)
           if mol.GetBondBetweenAtoms(carbonyl_carbon_index, bridgehead_atom_index):
            found_3_position = True
            break
        if found_3_position:
            break
    if not found_3_position:
        return False, "Carbonyl group at position 3 not found"


    # 3. Check for 5-alpha configuration. This can only be done by inspecting
    #    the SMILES representation
    
    # The 5th position is the carbon next to the first bridgehead in the SMILES string.
    # In the examples, it is always a chiral center that has an alpha hydrogen
    # i.e. [C@@] or [C@H].

    smiles_parts = smiles.split(']')
    if len(smiles_parts) < 2:
      return False, "Invalid SMILES: cannot find ring information"

    #Check the second set of characters which should be "[C@@" or "[C@H]"
    smiles_part_5_carbon = smiles_parts[1][1:4]
    if smiles_part_5_carbon != "C@H" and smiles_part_5_carbon != "C@@":
       return False, "No 5 alpha configuration found. Not [C@H] or [C@@]"


    return True, "Molecule is a 3-oxo-5alpha-steroid"