"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 is characterized by a sphingosine (d18:1) backbone, a phosphocholine
    head group, and a fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: True if molecule is a sphingomyelin d18:1, False otherwise, plus the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Sphingosine backbone pattern (C-C=C-C(OH)-C(NH)-C
    # This accounts for the stereochemistry
    sphingosine_core_pattern = Chem.MolFromSmarts("[C@@H]([C@H](C=C)CO)N[C@H]")
    if not mol.HasSubstructMatch(sphingosine_core_pattern):
         return False, "Sphingosine core not found"

    # 2. Phosphocholine head group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # 3. Fatty acyl chain (amide bond) pattern (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
         return False, "Fatty acyl chain not found"

    # 4.  Check for sphingosine d18:1  (18 carbons, one double bond in the sphingosine backbone)
    sphingosine_backbone_pattern = Chem.MolFromSmarts("CC[C@H]([C@H](C=C)CO)N[C@H]")
    matches = mol.GetSubstructMatches(sphingosine_backbone_pattern)
    if not matches:
        return False, "Sphingosine d18:1 backbone not found"

    # Locate the sphingosine carbon chain
    match = matches[0] # The pattern matches only once
    backbone_atoms = [match[0], match[1], match[2], match[3], match[4], match[5]]
    # Find an additional carbon attached to the backbone to count all the carbons
    found = False
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[0]).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms and neighbor.GetSymbol() == 'C':
            backbone_atoms.append(neighbor.GetIdx())
            found = True
            break
    if not found:
        return False, "Error identifying complete backbone chain."
    
    # Find all the carbons of the backbone
    last_carbon_idx = backbone_atoms[-1]
    while True:
      found = False
      for neighbor in mol.GetAtomWithIdx(last_carbon_idx).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms and neighbor.GetSymbol() == 'C':
          backbone_atoms.append(neighbor.GetIdx())
          last_carbon_idx = neighbor.GetIdx()
          found = True
          break
      if not found:
         break

    c_count = sum(1 for atom_idx in backbone_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Incorrect backbone carbon count. Should be 18, got {c_count}"

    # Count double bonds within the sphingosine
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            if bond.GetBeginAtomIdx() in backbone_atoms and bond.GetEndAtomIdx() in backbone_atoms:
                double_bond_count+=1
    if double_bond_count != 1:
        return False, f"Incorrect number of double bonds in backbone, should be 1, got {double_bond_count}"


    return True, "Meets all criteria for sphingomyelin d18:1"