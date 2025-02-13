"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is defined as a biomacromolecule consisting of >10 monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for repeating monosaccharide units (e.g. glucose, fructose)
    mon_smarts = ["OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", # Glucose
                  "O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)[C@@H]1O"] # Fructose
    mon_patterns = [Chem.MolFromSmarts(sma) for sma in mon_smarts]
    matches = [mol.GetSubstructMatches(pat) for pat in mon_patterns]
    if not any(match for match_list in matches for match in match_list):
        return False, "No monosaccharide units found"

    # Look for glycosidic linkages between monosaccharide units
    gly_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@@H]2O)CO)O[C@@H]1CO")
    gly_matches = mol.GetSubstructMatches(gly_pattern)
    if not gly_matches:
        return False, "No glycosidic linkages found"

    # Count number of monosaccharide units
    mon_count = len(matches[0])
    for match_list in matches[1:]:
        mon_count += len(match_list)
    
    if mon_count > 10:
        return True, f"Contains {mon_count} monosaccharide units linked glycosidically"
    else:
        return False, f"Only {mon_count} monosaccharide units found, need >10 for polysaccharide"