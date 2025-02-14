"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules arising from the oxidation of C20 EFAs.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for 20 Carbons in the molecule
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 20:
        return False, "Less than 20 carbons, not an icosanoid."

    # Check for C20 carbon backbone using substructure matching
    c20_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(c20_chain_pattern):
          return False, "No C20 backbone found."

    # 2. Check for Carboxylic Acid or Ester Group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")  # Corrected SMARTS for carboxylate
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxylic acid or ester group found."

    # 3. Check for oxygen atoms attached to C20 chain or in a ring system.
    oxygen_on_c20_pattern = Chem.MolFromSmarts("[CX4,CX3]~[OX2]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_on_c20_pattern)
    if len(oxygen_matches) == 0:
        # check if there is a ring with oxygens
        ring_with_oxygen = Chem.MolFromSmarts("[C]1[C]([O])[C][C][C]1")
        ring_oxygen_matches = mol.GetSubstructMatches(ring_with_oxygen)
        if len(ring_oxygen_matches) == 0:
          return False, "No oxygenation on C20 backbone or ring system"
    

    #4. Check for at least 2 double bonds - a strong indicator of Icosanoids, especially when conjugated
    num_double_bonds = 0
    for bond in mol.GetBonds():
          if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                num_double_bonds += 1
    if num_double_bonds < 2:
      return False, "Not enough double bonds for an Icosanoid"

    # 5. Check for conjugated system near oxygen/carbonyl/hydroxyl groups (more specific checks)

    # Check for prostane ring with oxygenation.
    prostane_pattern = Chem.MolFromSmarts("[C]1[C]([C])[C]([O])[C]([C])C1")
    prostane_ring_match = mol.GetSubstructMatches(prostane_pattern)
    if prostane_ring_match:
      # Check for prostane ring oxygenation:
      hydroxyl_on_ring = Chem.MolFromSmarts("[C]1[C]([O])[C][C][C]1")
      ketone_on_ring = Chem.MolFromSmarts("[C]1[C](=[O])[C][C][C]1")
      if mol.HasSubstructMatch(hydroxyl_on_ring) or mol.HasSubstructMatch(ketone_on_ring):
          return True, "Matches Icosanoid criteria (C20 backbone with prostane like ring with required oxygenation)"

    # Check for conjugated system and hydroxyl groups.
    conjugated_system_with_oxygen = Chem.MolFromSmarts("[C]=[C][C]=[C][C]([O])[C]")
    if mol.HasSubstructMatch(conjugated_system_with_oxygen):
      return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and oxygenation)"
    
    # Check for conjugated system, carbonyl and hydroxyl groups
    conjugated_system_with_carbonyl_and_hydroxyl = Chem.MolFromSmarts("[C]=[C][C]=[C][C](=O)[C](O)[C]")
    if mol.HasSubstructMatch(conjugated_system_with_carbonyl_and_hydroxyl):
      return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and carbonyl/hydroxyl group)"
      
    # Check for hydroperoxide on the C20 chain
    hydroperoxide_pattern = Chem.MolFromSmarts("[CX4,CX3]~[O]-[O]")
    if mol.HasSubstructMatch(hydroperoxide_pattern):
          return True, "Matches Icosanoid criteria (C20 backbone with hydroperoxide)"

    # Check for epoxides on the C20 chain
    epoxide_pattern = Chem.MolFromSmarts("[CX4,CX3]~C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
      return True, "Matches Icosanoid criteria (C20 backbone with epoxide)"
      
    return False, "Does not match Icosanoid criteria"