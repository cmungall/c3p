"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has a 3D conformation
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)  # Add hydrogens for conformer generation
        AllChem.EmbedMolecule(mol,randomSeed=42)
        if mol.GetNumConformers() == 0:
            return False, "Could not generate a 3D conformation for the molecule."


    # 1. Check for steroid core
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]([C]([C]3[C]([C]([C]1[C])([C]2)C)[C]4[C]([C](C3)([C])C)CC4)[H])([H])[H]")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Not a steroid core structure"

    # 2. Check for a hydroxy group at position 3 using SMARTS (ring position number)
    hydroxy_at_3_pattern = Chem.MolFromSmarts("[C;R][C;R]([O])([H])[C;R]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_at_3_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxyl group found at position 3"

    #Get the positions of the atoms in the substructure matching the steroid_core pattern:
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if len(core_match) != 20: #The pattern has 20 atoms.
        return False, "Internal error: wrong number of atoms in core"
    
    #Get the positions of the atoms in the substructure matching the hydroxy_at_3_pattern:
    hydroxy_match = mol.GetSubstructMatch(hydroxy_at_3_pattern)
    if len(hydroxy_match) != 3: #The pattern has 3 atoms.
        return False, "Internal error: wrong number of atoms in hydroxy_match"
    
    #3. Check for alpha-orientation of the hydroxyl group:
    conf = mol.GetConformer()
    
    #Indexes of atoms:
    pos_3 = hydroxy_match[1]
    pos_10 = core_match[12]
    pos_13 = core_match[15]
    pos_3_h = hydroxy_match[2]

    #Coordinates of atoms:
    p3_xyz = conf.GetAtomPosition(pos_3)
    p10_xyz = conf.GetAtomPosition(pos_10)
    p13_xyz = conf.GetAtomPosition(pos_13)
    p3h_xyz = conf.GetAtomPosition(pos_3_h)
    
    #Compute the centroid of positions 10 and 13
    p1013_center = (np.array(p10_xyz) + np.array(p13_xyz)) / 2.0
    
    #Direction vector from the centroid of 10/13 to H at pos 3
    dir_3h = np.array(p3h_xyz) - p1013_center

    #Direction vector from center 10/13 to O at pos 3.
    dir_3o = np.array(p3_xyz) - p1013_center

    #Check if the hydroxy group at position 3 is in the alpha configuration by checking the orientation of the H at position 3.
    dot_product = np.dot(dir_3h,dir_3o)
    if dot_product > 0: #If the dot product is positive, the H at pos 3 points in the same direction as the hydroxy group (not alpha)
        return False, "Hydroxyl group at position 3 is not in alpha-configuration"
        
    return True, "Molecule is a 3alpha-hydroxy steroid"