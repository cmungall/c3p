"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for flavan-3-ol unit (catechin/epicatechin core)
    # This pattern captures the chroman-3-ol core with a phenyl group at position 2
    flavan3ol_smarts = '[#6;R]1([#6])[#6;R][#8;R][#6;R]2[#6;R]1[#6]=[#6][#6]=[#6][#6]=2'
    flavan3ol = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol is None:
        return False, "Error in flavan-3-ol SMARTS pattern"
    
    # Find flavan-3-ol units in the molecule
    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol)
    num_flavan_units = len(flavan3ol_matches)
    
    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan-3-ol units, need at least 2"
    
    # Check for interflavan linkages between flavan-3-ol units
    # For each pair of flavan units, check if there is a bond between them
    bonds_between_units = 0
    for i in range(len(flavan3ol_matches)):
        for j in range(i+1, len(flavan3ol_matches)):
            unit_i_atoms = set(flavan3ol_matches[i])
            unit_j_atoms = set(flavan3ol_matches[j])
            # Check for bonds between atoms of unit_i and unit_j
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if (a1 in unit_i_atoms and a2 in unit_j_atoms) or (a2 in unit_i_atoms and a1 in unit_j_atoms):
                    bonds_between_units += 1
                    break  # One bond is enough to confirm linkage
    
    if bonds_between_units == 0:
        return False, "No interflavan bonds found between flavan-3-ol units"
    
    return True, "Contains two or more flavan-3-ol units connected via interflavan bonds"