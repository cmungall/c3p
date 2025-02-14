"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:33713 tannin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is defined as 'Any of a group of astringent polyphenolic vegetable principles or compounds,
    chiefly complex glucosides of catechol and pyrogallol.'
    This function checks for structural features characteristic of hydrolyzable and condensed tannins, 
    including galloyl groups attached to sugars and flavan-3-ol units linked via C-C or C-O bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for galloyl group attached via ester linkage to sugar (hydrolyzable tannins)
    galloyl_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1OC(=O)c2cc(O)c(O)c(O)c2')
    if mol.HasSubstructMatch(galloyl_sugar_pattern):
        return True, "Molecule is a gallotannin (hydrolyzable tannin with galloyl groups on sugar)"

    # Define pattern for hexahydroxydiphenoyl (HHDP) group attached to sugar (ellagitannins)
    hhdp_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O[C@H]2OC(=O)c3cc(O)c(O)c4oc(=O)c(O)cc4c3O[C@@H]2O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O')
    if mol.HasSubstructMatch(hhdp_sugar_pattern):
        return True, "Molecule is an ellagitannin (hydrolyzable tannin with HHDP group on sugar)"

    # Define pattern for flavan-3-ol units (catechin, epicatechin)
    flavan3ol_pattern = Chem.MolFromSmarts('OC[C@H]1Oc2cc(O)ccc2[C@@H]([C@H]1O)c1cc(O)cc(O)c1')
    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol_pattern)

    # Define patterns for interflavan linkages in condensed tannins (proanthocyanidins)
    # C4-C8 linkage
    c4_c8_linkage_pattern = Chem.MolFromSmarts('Oc1ccc(-c2c(O)cc3O[C@@H]4OC[C@H](O)[C@@H]4Oc3c2O)cc1')
    # C4-C6 linkage
    c4_c6_linkage_pattern = Chem.MolFromSmarts('Oc1cc(-c2c(O)cc3O[C@@H]4OC[C@H](O)[C@@H]4Oc3c2O)ccc1O')

    if len(flavan3ol_matches) >= 2 and (mol.HasSubstructMatch(c4_c8_linkage_pattern) or mol.HasSubstructMatch(c4_c6_linkage_pattern)):
        return True, "Molecule is a condensed tannin (proanthocyanidin) with flavan-3-ol units linked via interflavan bonds"

    # Check for phlorotannins (marine algae tannins composed of phloroglucinol units)
    phlorotannin_pattern = Chem.MolFromSmarts('Oc1cc(O)cc(O)c1')
    phlorotannin_matches = mol.GetSubstructMatches(phlorotannin_pattern)
    if len(phlorotannin_matches) >= 3:
        return True, "Molecule is a phlorotannin composed of multiple phloroglucinol units"

    # Check for complex tannins (combination of hydrolyzable and condensed tannin units)
    if (mol.HasSubstructMatch(galloyl_sugar_pattern) or mol.HasSubstructMatch(hhdp_sugar_pattern)) and len(flavan3ol_matches) >= 1:
        return True, "Molecule is a complex tannin containing both hydrolyzable and condensed tannin units"

    # Additional check for catechin/epicatechin gallate esters
    catechin_gallate_pattern = Chem.MolFromSmarts('OC(=O)c1cc(O)cc(O)c1')
    if len(flavan3ol_matches) >= 1 and mol.HasSubstructMatch(catechin_gallate_pattern):
        return True, "Molecule is a catechin or epicatechin gallate ester"

    # Check for high polyphenolic content (heuristic)
    num_phenol_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1cc(O)ccc1')))
    if num_phenol_groups >= 6 and rdMolDescriptors.CalcExactMolWt(mol) > 500:
        return True, f"Molecule has high polyphenolic content ({num_phenol_groups} phenol groups)"

    # If none of the above, molecule is likely not a tannin
    return False, "Molecule does not match structural features characteristic of tannins"