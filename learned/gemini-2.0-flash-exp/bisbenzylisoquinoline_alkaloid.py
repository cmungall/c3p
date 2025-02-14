"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid is characterized by two benzylisoquinoline units linked by ether bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
         str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for a benzylisoquinoline unit.
    # This pattern covers common substitution patterns around the isoquinoline and benzyl parts.
    # The aromatic ring can have X for C or N and includes substitutions. The rest is fixed.
    benzylisoquinoline_pattern = Chem.MolFromSmarts("[CX3,CX4,NX3,NX4]1=C2C=CC=CC=2C([CH2]C3=CC=CC=C3)=C1")
    if benzylisoquinoline_pattern is None:
        return None, "Invalid SMARTS Pattern"

    # Find all matches of the benzylisoquinoline pattern
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
         return False, f"Found {len(matches)} benzylisoquinoline units, need at least 2"

    # Check for ether linkages (C-O-C) connecting the two benzylisoquinoline units
    # the ether must connect the two rings
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    #  check if there is at least one ether
    if not ether_matches:
        return False, "No ether bridge connecting the benzylisoquinoline units"

    
    # Count the number of ether linkages
    ether_count = 0
    for match in ether_matches:
        atom1 = mol.GetAtomWithIdx(match[0])
        atom2 = mol.GetAtomWithIdx(match[1])
        atom3 = mol.GetAtomWithIdx(match[2])
    
        # Check if the neighboring atoms of the ether are part of a bis-benzylisoquinoline system
        is_connected_to_bzisoquinolines = False

        for bzisoq_match in matches:
            for bzisoq_idx in bzisoq_match:
               bzisoq_atom = mol.GetAtomWithIdx(bzisoq_idx)
               if atom1.GetIdx() in [neighbor.GetIdx() for neighbor in bzisoq_atom.GetNeighbors()]:
                   
                    for bzisoq_match2 in matches:
                        if bzisoq_match == bzisoq_match2:
                            continue
                        for bzisoq_idx2 in bzisoq_match2:
                             bzisoq_atom2 = mol.GetAtomWithIdx(bzisoq_idx2)
                             if atom3.GetIdx() in [neighbor.GetIdx() for neighbor in bzisoq_atom2.GetNeighbors()]:
                                 is_connected_to_bzisoquinolines = True
                                 break
                        if is_connected_to_bzisoquinolines:
                            break

            if is_connected_to_bzisoquinolines:
                break
        if is_connected_to_bzisoquinolines:
            ether_count+=1
        

    if ether_count < 1:
        return False, "The ether bridge does not connect the two benzylisoquinoline units"

    return True, "Contains two benzylisoquinoline units connected by at least one ether bridge"