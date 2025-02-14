"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: Glycolipid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as any molecule consisting of a glycosidic linkage
    between a carbohydrate moiety (usually a mono-, di-, or trisaccharide)
    and a lipid moiety (such as fatty acids or glycerolipids). In some cases,
    the glycerol backbone may be absent, and the sugar part may be acylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for carbohydrates (pyranose rings)
    sugar_pattern = Chem.MolFromSmarts("OC1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1")  # Generic hexose
    if sugar_pattern is None:
        return None, "Error in sugar SMARTS pattern"

    # Search for carbohydrate moiety
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No carbohydrate moiety found"

    # Define SMARTS pattern for long aliphatic chain (lipid moiety)
    lipid_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)")  # Generic long chain
    if lipid_pattern is None:
        return None, "Error in lipid SMARTS pattern"

    # Search for lipid moiety
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if len(lipid_matches) == 0:
        return False, "No lipid moiety found"

    # Check for glycosidic linkage between sugar and lipid
    # Find bonds connecting sugar and lipid moieties
    sugar_atoms = set()
    for match in sugar_matches:
        sugar_atoms.update(match)
    lipid_atoms = set()
    for match in lipid_matches:
        lipid_atoms.update(match)

    linked = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if (begin_atom in sugar_atoms and end_atom in lipid_atoms) or \
           (begin_atom in lipid_atoms and end_atom in sugar_atoms):
            # Check if the bond is an ether linkage (glycosidic bond)
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                linked = True
                break

    if linked:
        return True, "Contains carbohydrate and lipid moieties linked via glycosidic bond"
    else:
        return False, "No glycosidic linkage between carbohydrate and lipid moieties found"