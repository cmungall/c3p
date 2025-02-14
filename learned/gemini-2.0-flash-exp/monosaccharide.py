"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms and no glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Less than 3 carbon atoms"

    # Check for at least one carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group"
    

    # Check for glycosidic linkages (ether bond between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[CX4]~[OX2]~[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    for match in glycosidic_matches:
        oxygen_atom = mol.GetAtomWithIdx(match[1])
        if oxygen_atom.IsInRing():
            
            neighbor_carbon_1 = mol.GetAtomWithIdx(match[0])
            neighbor_carbon_2 = mol.GetAtomWithIdx(match[2])
            if neighbor_carbon_1.IsInRing() and neighbor_carbon_2.IsInRing() :
               return False, "Possible glycosidic linkage"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2][H]")  # -OH
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))


    return True, "Meets criteria for a monosaccharide"