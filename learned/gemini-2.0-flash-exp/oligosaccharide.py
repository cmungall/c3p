"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound formed by a few monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get number of rings
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()

    if n_rings < 2:
        return False, "Molecule has less than 2 rings, not an oligosaccharide"
    
    # Check for glycosidic linkages (C-O-C between rings)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[CX4][OX2]([CX4])")
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
      return False, "Molecule does not contain any glycosidic linkages."
    
    # Check that the molecule consists mainly of 5- and 6-membered rings with at least one oxygen
    #in the ring.
    sugar_ring_pattern1 = Chem.MolFromSmarts("[C;R1,R2]1[C;R1,R2][C;R1,R2][C;R1,R2][O;R1,R2]1") #5 member ring
    sugar_ring_pattern2 = Chem.MolFromSmarts("[C;R1,R2]1[C;R1,R2][C;R1,R2][C;R1,R2][C;R1,R2][O;R1,R2]1") #6 member ring
    
    if not mol.HasSubstructMatch(sugar_ring_pattern1) and not mol.HasSubstructMatch(sugar_ring_pattern2):
         return False, "Molecule does not contain sugar rings (5 or 6 membered cycles with an oxygen in the ring)."

    # Check all ring atoms are either C or O
    for atom in mol.GetAtoms():
      if atom.IsInRing() and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 8:
         return False, "Rings must contain only carbon and oxygen atoms."
    
    #Check the link between rings
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    for match in glycosidic_matches:
      start_atom = mol.GetAtomWithIdx(match[0])
      mid_atom = mol.GetAtomWithIdx(match[1])
      end_atom = mol.GetAtomWithIdx(match[2])
      if not (mid_atom.GetAtomicNum()==8 and start_atom.IsInRing() and end_atom.IsInRing()):
         return False, "Linkages between the cycles must be ether links between two rings."

    return True, "Molecule contains multiple sugar rings linked by glycosidic bonds."