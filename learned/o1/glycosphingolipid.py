"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide,
    where the carbohydrate is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ceramide pattern: sphingoid base with an amide linkage at C-2
    ceramide_smarts = "[C;H1,C;H2;H3][C@H](O)[C@H](NC(=O)[C;H1,C;H2;H3])[CH2][CH](O)[C;H1,C;H2;H3]"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if ceramide_pattern is None:
        return False, "Error in defining ceramide pattern"

    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide backbone found"

    # Define glycosidic linkage pattern: sugar attached via oxygen
    glycosidic_linkage_smarts = "[C@H]1([O,C])[O,C][C@H]([O,C])[C@@H]([O,C])[C@H]([O,C])[C@H]1O[*]"  # Sugar ring attached via O
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_linkage_smarts)
    if glycosidic_pattern is None:
        return False, "Error in defining glycosidic linkage pattern"

    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkage to sugar found"

    # Verify that the glycosidic linkage is connected to O-1 of the ceramide
    # Find the attachment point of the sugar
    attachment_found = False
    for match in ceramide_matches:
        ceramide_oxygen_idx = match[1]  # Index of the oxygen at C-1 in ceramide pattern
        for glyco_match in glycosidic_matches:
            sugar_attachment_idx = glyco_match[0]  # Index where sugar attaches
            bond = mol.GetBondBetweenAtoms(ceramide_oxygen_idx, sugar_attachment_idx)
            if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                attachment_found = True
                break
        if attachment_found:
            break

    if not attachment_found:
        return False, "Sugar moiety not attached via glycosidic linkage to O-1 of ceramide"

    return True, "Molecule is a glycosphingolipid with ceramide backbone and glycosidic linkage to sugar at O-1"