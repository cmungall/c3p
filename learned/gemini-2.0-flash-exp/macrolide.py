"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a macrocyclic lactone ring (12+ atoms) derived from polyketides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find lactone group
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)

    if not lactone_matches:
        return False, "No lactone group found"

    # Check for rings and if any ring contains a lactone group and has size 12+
    ring_info = mol.GetRingInfo()
    for ring_atom_ids in ring_info.AtomRings():
        if len(ring_atom_ids) >= 12:
            for lactone_match in lactone_matches:
                lactone_carbon_index = lactone_match[0]
                if lactone_carbon_index in ring_atom_ids:
                    # Count carbons in ring
                    carbon_count = sum(1 for atom_index in ring_atom_ids if mol.GetAtomWithIdx(atom_index).GetAtomicNum() == 6)
                    if carbon_count >= 10:
                        return True, f"Macrocyclic lactone ring (size: {len(ring_atom_ids)}, carbon count: {carbon_count}) found."

    return False, "No macrocyclic lactone ring (12+ atoms) containing a lactone group found with minimum 10 carbon atoms."