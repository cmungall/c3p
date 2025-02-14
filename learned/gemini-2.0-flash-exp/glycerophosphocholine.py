"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: CHEBI:17797 glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine consists of a glycerol backbone, a phosphate group, a phosphocholine headgroup,
    and typically two fatty acid chains attached via ester or ether bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol-phosphate fragment, with a free connection for the fatty acid chains.
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]([OX2])([CHX4]([OX2])[CH2X4][OX2][PX4](=[OX1])[OX1])")

    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol phosphate backbone detected"

    # Define the phosphocholine head group pattern connected to the phosphate
    phosphocholine_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
          return False, "No phosphocholine group detected"
    
    # Check for at least one acyl chain:
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 1:
       # Check ether-linked chains
       ether_pattern = Chem.MolFromSmarts("[OX2][CX4,CX3]")
       ether_matches = mol.GetSubstructMatches(ether_pattern)
       if len(ether_matches) < 1:
            return False, "No fatty acid or ether-linked chain found"
       
    # Check for at least one fatty acid chain of adequate length (at least 4 carbons long)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "No fatty acid chain found"

    # Count carbons, oxygens, phosphorus and nitrogen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if p_count != 1:
        return False, f"Incorrect number of phosphates, found {p_count}"
    if n_count != 1:
        return False, f"Incorrect number of nitrogens, found {n_count}"
    
    if c_count < 5:
        return False, "Too few carbons for glycerophosphocholine"
    
    if o_count < 4:
      return False, "Too few oxygens for glycerophosphocholine"
    
    return True, "Glycerophosphocholine structure detected"