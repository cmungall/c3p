"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is a sphingoid base with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Sphingosine Backbone Pattern
    # Focus on core with a double bond on one side, and a -CH2OH on the other, and a -OH next to it
    # the double bond may be on either side of the central part
    sphingosine_pattern1 = Chem.MolFromSmarts("[CX4]=[CX3]-[CX4]([OH])-[CX4]([NH])-[CX4]([CH2X4][OH])") #double bond, C-OH-C-N-C-CH2OH
    sphingosine_pattern2 = Chem.MolFromSmarts("[CX4]([CH2X4][OH])-[CX4]([NH])-[CX4]([OH])-[CX3]=[CX4]") #CH2OH-C-N-C-OH-C double bond
    if not (mol.HasSubstructMatch(sphingosine_pattern1) or mol.HasSubstructMatch(sphingosine_pattern2)):
        return False, "No sphingosine backbone found"
        
    # 2. Amide Linkage Pattern
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if not amide_matches:
        return False, "No amide linkage found"

    # Check that the amide is linked to the sphingosine
    found_amide_connected_to_sphingosine = False
    for amide_match in amide_matches:
        for atom_index in amide_match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == 'N':
                for neighbor in atom.GetNeighbors():
                     if neighbor.Match(Chem.MolFromSmarts("[CX4]([OH])-[CX4]([NH])-[CX4]([CH2X4][OH])")) or \
                           neighbor.Match(Chem.MolFromSmarts("[CX4]([CH2X4][OH])-[CX4]([NH])-[CX4]([OH])")) or \
                           neighbor.Match(Chem.MolFromSmarts("[CX4]=[CX3]-[CX4]([OH])-[CX4]([NH])")) or \
                           neighbor.Match(Chem.MolFromSmarts("[CX4]([OH])-[CX3]=[CX4]")):
                        found_amide_connected_to_sphingosine = True
                        break
            if found_amide_connected_to_sphingosine:
              break
        if found_amide_connected_to_sphingosine:
          break
    
    if not found_amide_connected_to_sphingosine:
        return False, "Amide not connected to the sphingosine backbone"


    # 3. Fatty Acid Chain Check (using rotatable bonds)
    # We count rotatable bonds connected to carbonyl
    for match in amide_matches:
      carbonyl_c_index = -1
      for atom_index in match:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 3 and any (neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors()):
            carbonyl_c_index = atom_index
            break
      if carbonyl_c_index == -1:
          continue
      submol = Chem.Mol(mol)
      for atom in submol.GetAtoms():
            if atom.GetIdx() != carbonyl_c_index:
                 submol.GetAtomWithIdx(atom.GetIdx()).SetAtomMapNum(0)

      submol = Chem.DeleteSubstructs(submol,Chem.MolFromSmarts("[#0]"))
      
      num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(submol)
      if num_rotatable_bonds < 12 or num_rotatable_bonds > 25:
          return False, f"Fatty acid chain length not within range (14-26 C). Found {num_rotatable_bonds}"


    return True, "Contains a sphingosine backbone with an amide-linked fatty acid"