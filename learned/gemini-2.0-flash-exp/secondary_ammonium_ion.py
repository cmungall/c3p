"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is a nitrogen atom with a positive charge and two directly attached
    carbon containing groups via single bonds, and *none* of these carbons are part of any ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    has_secondary_ammonium = False
    
    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1: # Check for nitrogen and positive charge
            carbon_neighbors = 0
            neighbor_carbons = []
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == rdchem.BondType.SINGLE: # Check if neighbour is a carbon via single bond
                    carbon_neighbors += 1
                    neighbor_carbons.append(neighbor)
            if carbon_neighbors == 2:
              all_carbons_outside_ring = True
              for carbon_neighbor in neighbor_carbons:
                if carbon_neighbor.IsInRing():
                  all_carbons_outside_ring = False
                  break
              if all_carbons_outside_ring:
                has_secondary_ammonium = True
                break
            
    if has_secondary_ammonium:
         return True, "Contains at least one positively charged nitrogen with two directly attached carbon groups via single bonds, and *none* of these carbons are part of any ring."
    else:
        return False, "Does not contain a positively charged nitrogen with two directly attached carbon groups via single bonds, and *none* of these carbons are part of any ring."