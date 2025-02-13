"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is of the acrovestone class.
    Acrovestone is characterized as a polyphenol, typically isolated from a plant species, with antioxidant activities. 
    Common structural features include isoflavonoid cores and glycosidic attachments.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as acrovestone, False otherwise
        str: Reason for classification
    """
    
    # Try to parse the SMILES string, check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for polyphenolic structure by searching for aromatics with hydroxyl groups
    # SMARTS pattern: a set of aromatic rings with hydroxyl groups attached at different positions 
    polyphenol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    if not mol.HasSubstructMatch(polyphenol_pattern):
        return False, "No polyphenolic structure detected"

    # 2. Check for glycosylation (presence of glycosidic bonds)
    # A common pattern for a sugar moiety features -O-CO- group attached, indicating a glucoside connection
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "Glycosylation (such as glucoside) not detected"

    # Because glycosides are prominent, ensure glycosidic linkage is present
    # Iterate over glycoside patterns if necessary

    # Further verification could involve counting the number of aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, "Too few aromatic rings for a polyphenolic structure"

    # Additional attributes specific to the plant origin may need mass spectral analysis

    return True, "Molecule is consistent with structural features of acrovestone"