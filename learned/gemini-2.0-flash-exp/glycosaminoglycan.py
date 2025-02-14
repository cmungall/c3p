"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for presence of aminomonosaccharide rings using a more general SMARTS pattern.
    # This pattern looks for a ring with a nitrogen atom directly bonded to a carbon of the ring
    # The nitrogen can be part of an amide (acetylated) or a free amine
    amino_sugar_pattern = Chem.MolFromSmarts("[NX3;!H0][CX4;R]")
    if amino_sugar_pattern is None:
         return False, "Invalid SMARTS pattern for aminosugar."
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)

    if not amino_sugar_matches:
        return False, "No aminomonosaccharide unit detected."

    # 2. Check for glycosidic linkages between two ring carbons.
    # This pattern looks for a carbon of ring linked to an oxygen, which in turn is linked to another carbon
    glycosidic_link_pattern = Chem.MolFromSmarts("[CX4;R]~[OX2]~[CX4;R]")
    if glycosidic_link_pattern is None:
          return False, "Invalid SMARTS pattern for glycosidic linkage."
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_link_pattern)

    if len(glycosidic_matches) == 0:
         return False, "No glycosidic linkage between rings detected."

    # 3. Minimal molecular weight check - GAGs are polysaccharides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a polysaccharide"

    # 4. Check for more than one ring
    ring_pattern = Chem.MolFromSmarts("[R]")
    if ring_pattern is None:
         return False, "Invalid SMARTS pattern for ring."
    ring_matches = mol.GetSubstructMatches(ring_pattern)
    if len(ring_matches) <= 1:
         return False, "Need at least two rings for a polysaccharide."


    return True, "Contains aminomonosaccharide units linked by glycosidic bonds."