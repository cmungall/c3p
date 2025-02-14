"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone with correct stereochemistry (sn-glycerol)
    # The middle carbon should have configuration R.
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H]([CH2X4])")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
       return False, "No glycerol backbone with correct chirality found"
    glycerol_c2 = glycerol_match[1]
    
    # 2. Check for phosphocholine head group
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)(O)(OCC[N+](C)(C)C)[O-]")
    phosphocholine_match = mol.GetSubstructMatch(phosphocholine_pattern)
    if not phosphocholine_match:
        return False, "No phosphocholine group found"
    
    # Verify phosphocholine is attached to glycerol
    phospho_c = mol.GetAtomWithIdx(phosphocholine_match[2])
    glycerol_c1 = glycerol_match[0]
    glycerol_c3 = glycerol_match[2]

    glycerol_c1_atom = mol.GetAtomWithIdx(glycerol_c1)
    glycerol_c3_atom = mol.GetAtomWithIdx(glycerol_c3)


    phospho_attached_to_glycerol = False
    for neighbor in glycerol_c3_atom.GetNeighbors():
      if neighbor.GetIdx() == phospho_c.GetIdx():
            phospho_attached_to_glycerol = True
            break

    if not phospho_attached_to_glycerol:
       return False, "Phosphocholine not attached to glycerol"


    # 3. Check for ether linkage at position 1
    ether_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    ether_found = False
    for match in ether_matches:
        if match[0] == glycerol_c1:
            ether_found = True
            ether_c = match[2]
            break
    if not ether_found:
        return False, "No ether linkage at position 1 found"

    # 4. Check for ester linkage at position 2
    ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_found = False
    for match in ester_matches:
      if match[0] == glycerol_c2:
        ester_found=True
        ester_c = match[2]
        break
    if not ester_found:
      return False, "No ester linkage at position 2 found"


    # 5. Check for alkyl and acyl chains (must contain long carbon chains)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    alkyl_chain_matches = []
    for atom in mol.GetAtomWithIdx(ether_c).GetNeighbors():
      temp_mol = Chem.MolFromSmarts(f"[{atom.GetSymbol()}]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
      alkyl_matches = mol.GetSubstructMatches(temp_mol)
      if len(alkyl_matches) > 0:
        alkyl_chain_matches += alkyl_matches

    if len(alkyl_chain_matches) < 1:
        return False, "No long alkyl chain at position 1"


    acyl_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])~[CX4]~[CX4]~[CX4]~[CX4]")
    acyl_chain_matches = []
    for atom in mol.GetAtomWithIdx(ester_c).GetNeighbors():
      temp_mol = Chem.MolFromSmarts(f"[{atom.GetSymbol()}]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
      acyl_matches = mol.GetSubstructMatches(temp_mol)
      if len(acyl_matches) > 0:
        acyl_chain_matches += acyl_matches

    if len(acyl_chain_matches) < 1:
       return False, "No long acyl chain at position 2"

    # Check chain lengths:
    alkyl_rotatable_bonds = 0
    acyl_rotatable_bonds = 0
    
    
    alkyl_chain_atoms = []
    queue = [mol.GetAtomWithIdx(ether_c)]
    visited = {mol.GetAtomWithIdx(ether_c).GetIdx()}
    while len(queue) > 0:
      current_atom = queue.pop(0)
      alkyl_chain_atoms.append(current_atom)
      for neighbor in current_atom.GetNeighbors():
          if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum()==6:
            queue.append(neighbor)
            visited.add(neighbor.GetIdx())
    
    alkyl_rotatable_bonds = 0
    for atom in alkyl_chain_atoms:
      for bond in atom.GetBonds():
          if bond.GetBeginAtom().GetIdx() in visited and bond.GetEndAtom().GetIdx() in visited and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
              alkyl_rotatable_bonds += 1

    acyl_chain_atoms = []
    queue = [mol.GetAtomWithIdx(ester_c)]
    visited = {mol.GetAtomWithIdx(ester_c).GetIdx()}
    while len(queue) > 0:
      current_atom = queue.pop(0)
      acyl_chain_atoms.append(current_atom)
      for neighbor in current_atom.GetNeighbors():
          if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum()==6:
            queue.append(neighbor)
            visited.add(neighbor.GetIdx())

    acyl_rotatable_bonds = 0
    for atom in acyl_chain_atoms:
      for bond in atom.GetBonds():
          if bond.GetBeginAtom().GetIdx() in visited and bond.GetEndAtom().GetIdx() in visited and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
             acyl_rotatable_bonds += 1

    if alkyl_rotatable_bonds < 4 or acyl_rotatable_bonds < 4:
         return False, "Alkyl and/or Acyl chains too short."


    return True, "Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine criteria"