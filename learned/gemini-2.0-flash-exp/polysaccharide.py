"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a polymer of monosaccharides linked by glycosidic bonds,
    with more than 10 monosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Relaxed SMARTS pattern for a generic monosaccharide ring (any ring with at least one oxygen and at least 4 atoms).
    monosaccharide_ring_smarts = "[CX3,CX4;!H0;R]~[CX3,CX4;!H0;R]~[CX3,CX4;!H0;R]~[CX3,CX4;!H0;R]~[OX2;R]"
    monosaccharide_pattern = Chem.MolFromSmarts(monosaccharide_ring_smarts)
    if monosaccharide_pattern is None:
        return None, "Invalid SMARTS pattern, check syntax"

    # Find all matches of the monosaccharide ring pattern
    matches = mol.GetSubstructMatches(monosaccharide_pattern)
    num_rings = len(matches)

    # Count glycosidic bonds (C-O-C linkages) between rings. The carbons MUST be part of a ring
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4;!H0;R]-[OX2]-[CX4;!H0;R]")  #C-O-C linkage where Carbons are part of a ring
    if glycosidic_bond_pattern is None:
        return None, "Invalid SMARTS pattern for glycosidic bonds"
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic = len(glycosidic_matches)

    if num_rings == 0:
      return False, f"No monosaccharide rings found. {num_rings} rings, {num_glycosidic} bonds"

    if num_glycosidic == 0 and num_rings > 1 :
      return False, f"Multiple monosaccharide rings detected, but no glycosidic bonds. {num_rings} rings, {num_glycosidic} bonds"
    
    num_monosaccharide_units = num_glycosidic + 1
    if num_rings >= 10 or num_monosaccharide_units >=10:
      return True, f"Contains {num_rings} monosaccharide rings and at least {num_monosaccharide_units} glycosidic bonds. Likely a polysaccharide."
    
    # Check molecular weight - polysaccharides typically > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, f"Molecular weight too low for polysaccharide. {num_rings} rings, {num_glycosidic} bonds"
      
    # Heuristic approach:
    if num_rings > 3 and num_glycosidic > 1:
      return True, f"Multiple monosaccharide rings, linked by glycosidic bonds, suggesting an oligosaccharide/polysaccharide. {num_rings} rings, {num_glycosidic} bonds"

    # If all checks passed, but still below 10
    if num_rings > 2 and num_glycosidic > 1:
      return True, f"Multiple monosaccharide rings, possibly an oligosaccharide. Needs further analysis. {num_rings} rings, {num_glycosidic} bonds"

    return False, f"Fails polysaccharide criteria. {num_rings} rings, {num_glycosidic} bonds. Needs further analysis."