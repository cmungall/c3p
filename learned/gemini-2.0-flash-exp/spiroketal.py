"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if there are at least two rings in the molecule
    if mol.GetRingInfo().NumRings() < 2:
         return False, "Molecule does not have at least two rings."


    # 1. Identify ketal carbons: carbon bonded to two oxygens (sp3 hybridized), zero hydrogens and two carbons (within rings).
    ketal_pattern = Chem.MolFromSmarts("[CX4H0](-[OX2;R])(-[OX2;R])([CX4;R])([CX4;R])")
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)

    if not ketal_matches:
        return False, "No ketal carbon found."
    
    
    for match in ketal_matches:
        ketal_carbon_idx = match[0]
        ketal_carbon = mol.GetAtomWithIdx(ketal_carbon_idx)
        
        # Get oxygen atoms bonded to ketal carbon
        oxygen_indices = [neighbor.GetIdx() for neighbor in ketal_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 8]
        if len(oxygen_indices) != 2:
            continue # Should not happen

        oxy1_idx = oxygen_indices[0]
        oxy2_idx = oxygen_indices[1]

        # Verify that both oxygens are part of rings.
        oxy1 = mol.GetAtomWithIdx(oxy1_idx)
        oxy2 = mol.GetAtomWithIdx(oxy2_idx)
        if not oxy1.IsInRing() or not oxy2.IsInRing():
            continue

        # 2. Find the smallest ring containing each oxygen
        ring_info = mol.GetRingInfo()
        
        rings_oxy1 = [set(ring) for ring in ring_info.AtomRings() if oxy1_idx in ring]
        rings_oxy2 = [set(ring) for ring in ring_info.AtomRings() if oxy2_idx in ring]

        if not rings_oxy1 or not rings_oxy2:
             continue


        # Find smallest rings by size
        if rings_oxy1 and rings_oxy2:
            smallest_ring_oxy1 = min(rings_oxy1, key=len)
            smallest_ring_oxy2 = min(rings_oxy2, key=len)


            # 3. Verify if the two smallest rings intersect ONLY at the ketal carbon
            smallest_ring_oxy1.discard(ketal_carbon_idx)
            smallest_ring_oxy2.discard(ketal_carbon_idx)
            if not smallest_ring_oxy1.intersection(smallest_ring_oxy2):
                return True, "Spiroketal structure identified."
            
    return False, "No spiroketal structure detected."