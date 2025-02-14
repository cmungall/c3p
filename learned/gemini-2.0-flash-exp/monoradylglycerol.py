"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with a single acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 OH and 1 attachment)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Find the attachment points on the glycerol
    attachment_pattern = Chem.MolFromSmarts("[CH2X4;!H0]([OX2H1])[CHX4;!H0]([OX2H1])[CH2X4;!H0]")
    attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
       attachment_pattern = Chem.MolFromSmarts("[CH2X4;!H0]([OX2H1])[CHX4;!H0][CH2X4;!H0]([OX2H1])")
       attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
        attachment_pattern = Chem.MolFromSmarts("[CH2X4;!H0][CHX4;!H0]([OX2H1])[CH2X4;!H0]([OX2H1])")
        attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
        return False, "No attachment point found"

    # Count ester linkages (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Count ether linkages (-C-O-C-)
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Count vinyl ether linkages (-C=C-O-C-)
    vinylether_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OX2][CX4]")
    vinylether_matches = mol.GetSubstructMatches(vinylether_pattern)


    total_attachments = len(ester_matches) + len(ether_matches) + len(vinylether_matches)

    if total_attachments != 1:
        return False, f"Found {total_attachments} attachments, require exactly 1"

    # Count number of carbons and oxygens for consistency
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 3:
        return False, "Too few carbons"
    if o_count < 3:
         return False, "Too few oxygens"
    
    return True, "Contains a glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent"