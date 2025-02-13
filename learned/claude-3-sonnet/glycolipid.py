"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:16321 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a molecule with a carbohydrate part connected to a lipid part via a glycosidic linkage.

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

    # Look for carbohydrate patterns
    carbohydrate_patterns = [
        Chem.MolFromSmarts("[OX2][C@H][C@H](O)[C@H](O)[C@@H](O)C"),  # glucose
        Chem.MolFromSmarts("[OX2][C@H][C@H](O)[C@H](O)[C@H](O)C"),   # galactose
        Chem.MolFromSmarts("[OX2][C@H][C@H](O)[C@H](O)[C@H](O)C(O)"), # glucuronic acid
        # Add more patterns for other monosaccharides and oligosaccharides
    ]
    carbohydrate_matches = [mol.HasSubstructMatch(pattern) for pattern in carbohydrate_patterns]
    if not any(carbohydrate_matches):
        return False, "No carbohydrate part found"

    # Look for lipid part patterns
    lipid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCC"),  # Linear aliphatic chain
        Chem.MolFromSmarts("CC(C)CCCCC"),   # Branched aliphatic chain
        Chem.MolFromSmarts("C1CCCCCCCCC1"),  # Cyclic aliphatic chain
        # Add more patterns for other lipid moieties
    ]
    lipid_matches = [mol.HasSubstructMatch(pattern) for pattern in lipid_patterns]
    if not any(lipid_matches):
        return False, "No lipid part found"

    # Look for glycosidic linkage patterns
    glycosidic_linkage_patterns = [
        Chem.MolFromSmarts("[OX2][C@H][C@H](O)[C@H](O)[C@@H](O)CO"),  # Glycosidic oxygen linked to carbohydrate
        Chem.MolFromSmarts("[OX2][C@H][C@H](O)[C@H](O)[C@H](O)COC"),  # Glycosidic oxygen linked to lipid
        # Add more patterns for other glycosidic linkages
    ]
    glycosidic_linkage_matches = [mol.HasSubstructMatch(pattern) for pattern in glycosidic_linkage_patterns]
    if not any(glycosidic_linkage_matches):
        return False, "No glycosidic linkage found"

    # Look for specific glycolipid classes
    glycolipid_class_patterns = [
        Chem.MolFromSmarts("C[C@@H](O)[C@@H](O)COC[C@H](O)[C@H](O)CO"),  # Cerebroside
        Chem.MolFromSmarts("[O-]C(=O)C[C@@H](O)[C@H](O)COC[C@H](O)[C@H](O)CO"),  # Ganglioside
        # Add more patterns for other glycolipid classes
    ]
    glycolipid_class_matches = [mol.HasSubstructMatch(pattern) for pattern in glycolipid_class_patterns]
    if any(glycolipid_class_matches):
        return True, "Matches a known glycolipid class"

    # Check molecular weight and atom count ratios
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 2000:
        return False, "Molecular weight outside typical range for glycolipids"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 20 or o_count < 5:
        return False, "Insufficient carbon or oxygen atoms for a glycolipid"

    return True, "Contains a carbohydrate part connected to a lipid part via a glycosidic linkage"