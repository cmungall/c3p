"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a core trisaccharide, oligosaccharide side chains and
    3-hydroxytetradecanoic acid units, and also often phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for heptose (7-carbon sugar) units:
    # SMARTS pattern for a 7-carbon sugar, without checking stereochemistry
    heptose_pattern = Chem.MolFromSmarts("C1([CH]([CH]([CH]([CH]([CH](C1O)O)O)O)O)O)CO")
    heptose_matches = mol.GetSubstructMatches(heptose_pattern)

    # 2. Check for octulosonic acid (8-carbon sugar with keto group):
    # SMARTS pattern for octulosonic acid
    octulosonic_acid_pattern = Chem.MolFromSmarts("C(=O)C([CH]([CH]([CH]([CH]([CH](CO)O)O)O)O)O)O")
    octulosonic_acid_matches = mol.GetSubstructMatches(octulosonic_acid_pattern)


    # Check for at least one heptose and one octulosonic acid
    if len(heptose_matches) < 1 or len(octulosonic_acid_matches) < 1:
          return False, f"Not a lipopolysaccharide: Insufficient heptose or octulosonic acid units. found {len(heptose_matches)} heptose and {len(octulosonic_acid_matches)} octulosonic acid"


    # 3. Check for 3-hydroxytetradecanoic acid units (-O-C(=O)-CH(O)-C12 chain)
    hydroxytetradecanoic_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4](O)[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    hydroxytetradecanoic_matches = mol.GetSubstructMatches(hydroxytetradecanoic_pattern)
    if len(hydroxytetradecanoic_matches) < 1:
      return False, f"Not a lipopolysaccharide: Missing hydroxytetradecanoic acid units. found {len(hydroxytetradecanoic_matches)}"


    # 4. Check for Phosphate groups ( P(=O)(O)O or variations of it like P(=O)(O)[O-])
    phosphate_pattern1 = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_pattern2 = Chem.MolFromSmarts("P(=O)(O)[O-]")
    phosphate_matches1 = mol.GetSubstructMatches(phosphate_pattern1)
    phosphate_matches2 = mol.GetSubstructMatches(phosphate_pattern2)
    if (len(phosphate_matches1) + len(phosphate_matches2) ) < 1:
        return False, f"Not a lipopolysaccharide: Missing phosphate groups. found {len(phosphate_matches1)+len(phosphate_matches2)}"


    return True, "Lipopolysaccharide features found: heptose and octulosonic acid units, hydroxytetradecanoic acids and phosphate groups."