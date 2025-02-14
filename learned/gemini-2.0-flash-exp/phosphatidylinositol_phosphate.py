"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate has a glycerol backbone, a phosphate group,
    a myo-inositol ring, and fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone (C-C-C with an oxygen attached to each C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]") # Changed this from Triglyceride since the 3 carbons have an O attachement
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for a phosphate group attached to glycerol. Specifically look for -C-O-P(=O)(O)-O
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[PX4](=[OX1])(-[OX2])-[OX2]")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
          return False, "No phosphate attached to the glycerol"

    # 3. Check for myo-inositol ring (six-membered ring with 6 hydroxyls). Can be difficult with SMARTS directly since not all the OH groups are directly attached
    inositol_pattern = Chem.MolFromSmarts("C1C(C(C(C(C1O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # 4. Count phosphates attached to myo-inositol. This can range from 1 to 3 phosphates.
    # Find the inositol ring first (C1(C(C(C(C(C1O)O)O)O)O)O),
    # then find how many phosphorus atoms are attached to this ring
    inositol_match = mol.GetSubstructMatch(inositol_pattern)
    phosphates_on_inositol = 0
    if inositol_match:
        for atom_index in inositol_match:
            atom = mol.GetAtomWithIdx(atom_index)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 15: # Phosphorus
                  phosphates_on_inositol += 1
    if phosphates_on_inositol < 1: # At least one phosphate must be attached
        return False, f"Found {phosphates_on_inositol} phosphate groups on inositol, need at least 1"

    # 5. Check for ester groups (-O-C(=O)-) and fatty acid chains
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
      return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains. Use a looser definition than triglyceride
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
       return False, "Chains too short to be fatty acids"

    # Molecular weight - typical range is > 700
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
       return False, "Molecular weight too low for PIP"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 20:
       return False, "Too few carbons"
    if p_count < 1:
      return False, "Must have at least one phosphorus"
    if o_count < 10:
      return False, "Too few oxygen atoms"


    return True, "Contains a glycerol backbone, a phosphate, an inositol ring with at least one phosphate, and fatty acid chains"