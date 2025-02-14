"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone with fatty acid chains attached
    at positions sn-1 and sn-2 via ester or ether bonds, and a phosphocholine group at position sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Glycerol backbone with positions sn-1, sn-2, sn-3
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2][CH](O)[CH2][O]")
    if glycerol_pattern is None:
        return False, "Invalid glycerol backbone SMARTS pattern"

    # Phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("O[P](=O)([O-])OCC[N+](C)(C)C")
    if phosphocholine_pattern is None:
        return False, "Invalid phosphocholine SMARTS pattern"

    # Check for glycerol backbone
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Check for phosphocholine group
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for fatty acid chains attached via ester or ether bonds at sn-1 and sn-2
    # Ester or ether bonds can be represented by O-C pattern
    ester_ether_pattern = Chem.MolFromSmarts("O[C]")
    if ester_ether_pattern is None:
        return False, "Invalid ester/ether SMARTS pattern"

    fatty_acid_chains = 0
    # Check attachments at sn-1 and sn-2 positions in the glycerol backbone
    for match in glycerol_matches:
        sn1_O_idx = match[0]  # Oxygen at sn-1
        sn2_O_idx = match[2]  # Oxygen at sn-2

        # Check for ester or ether at sn-1
        sn1_O_atom = mol.GetAtomWithIdx(sn1_O_idx)
        sn1_attached = False
        for bond in sn1_O_atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(sn1_O_atom)
            if neighbor_atom.GetIdx() != match[1]:  # Exclude bond to glycerol backbone
                if mol.GetSubstructMatch(ester_ether_pattern, useChirality=False):
                    fatty_acid_chains += 1
                    sn1_attached = True
                    break

        # Check for ester or ether at sn-2
        sn2_O_atom = mol.GetAtomWithIdx(sn2_O_idx)
        sn2_attached = False
        for bond in sn2_O_atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(sn2_O_atom)
            if neighbor_atom.GetIdx() != match[1]:  # Exclude bond to glycerol backbone
                if mol.GetSubstructMatch(ester_ether_pattern, useChirality=False):
                    fatty_acid_chains += 1
                    sn2_attached = True
                    break

        if fatty_acid_chains >= 2:
            return True, "Molecule is a glycerophosphocholine"

    return False, "Does not match glycerophosphocholine structure"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'glycerophosphocholine',
        'definition': 'The glycerol phosphate ester of a phosphocholine. A nutrient with many different roles in human health.',
        'parents': []
    }
}