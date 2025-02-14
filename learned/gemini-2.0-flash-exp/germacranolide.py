"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are sesquiterpene lactones with a 10-membered ring and a lactone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the presence of a 10-membered ring
    has_10_ring = False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 10:
            has_10_ring = True
            break

    if not has_10_ring:
        return False, "No 10-membered ring found"
    
    # 2. Check for the presence of a lactone ring (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C(=O)O[C;R]") # Check for C=O-O within a ring
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if len(lactone_matches) != 1:
        return False, "Incorrect number of lactone rings found"

    # 3. Verify that lactone is connected to 10-member ring:
    lactone_match_atoms = set(lactone_matches[0])
    lactone_connected_to_10_ring = False

    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 10:
            ten_ring_atoms = set(ring)
            if len(lactone_match_atoms.intersection(ten_ring_atoms)) > 0:
                lactone_connected_to_10_ring = True
                break
    if not lactone_connected_to_10_ring:
      return False, "Lactone ring not connected to the 10-membered ring"
    

    # 4. Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 350:  # Sesquiterpenes are typically ~C15
        return False, f"Molecular weight ({mol_wt:.2f}) outside typical sesquiterpene range"

    # 5. Ring count
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Must have at least two rings"


    # If all checks pass, then classify as germacranolide
    return True, "Contains a 10-membered ring and a lactone connected to it, consistent with a germacranolide"