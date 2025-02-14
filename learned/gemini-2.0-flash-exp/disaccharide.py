"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharides linked by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of ring systems
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
       return False, f"Molecule has less than 2 rings"


    # Check for glycosidic bonds (C-O-C connecting two rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]O[CX4]")
    matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    ring_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]

    valid_glycosidic_bonds = 0
    for match in matches:
        
        atom1_idx = match[0]
        atom2_idx = match[2]
        if (atom1_idx in ring_atoms) and (atom2_idx in ring_atoms):
           valid_glycosidic_bonds += 1
        
    if valid_glycosidic_bonds < 1:
        return False, "No glycosidic bond connecting two rings found."

    # Verify a limited number of C and O atoms to avoid long chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 8 or c_count > 24:
      return False, "Carbon count outside expected range."
    if o_count < 4 or o_count > 12:
      return False, "Oxygen count outside expected range."


    return True, "Contains two monosaccharide units linked by a glycosidic bond"