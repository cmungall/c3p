"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached 
to an unspecified sugar via an anomeric diphosphate linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for uracil ring
    uracil_pattern = Chem.MolFromSmarts("O=c1cc[nH]c(=O)[nH]1")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil ring found"
    
    # Check for ribose attached to uracil (uridine)
    uridine_pattern = Chem.MolFromSmarts("O=c1cc[nH]c(=O)n1[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine moiety found"
        
    # Check for diphosphate group
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"
        
    # Check for sugar characteristics (multiple OH groups, ring oxygen)
    sugar_pattern = Chem.MolFromSmarts("[OH1,OH0][C@@H,C@H,CH1,CH2]1[O][C@@H,C@H,CH1,CH2][C@@H,C@H,CH1,CH2][C@@H,C@H,CH1,CH2][C@@H,C@H,CH1,CH2]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
        
    # Check that sugar is connected via anomeric carbon to diphosphate
    # This is a simplified pattern that looks for C-O-P linkage in a ring system
    anomeric_linkage = Chem.MolFromSmarts("[C;R1]([O;R1,r5,r6])OP(=O)(O)O")
    if not mol.HasSubstructMatch(anomeric_linkage):
        return False, "Sugar not connected via anomeric carbon to diphosphate"
        
    # Count phosphorus atoms - should be exactly 2 for UDP
    p_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[P]")))
    if p_count != 2:
        return False, f"Found {p_count} phosphorus atoms, need exactly 2"
        
    # Verify presence of key functional groups typical in UDP sugars
    if not (mol.HasSubstructMatch(Chem.MolFromSmarts("O[C@@H,C@H,CH1]1O[C@@H,C@H,CH1]")) and  # Ring oxygen with adjacent carbons
            mol.HasSubstructMatch(Chem.MolFromSmarts("COP(=O)(O)O")) and  # Phosphate connection
            mol.HasSubstructMatch(Chem.MolFromSmarts("[C@@H,C@H,CH1]O[C@@H,C@H,CH1]"))): # C-O-C linkages
        return False, "Missing key structural features of UDP-sugar"
        
    return True, "Contains UDP moiety connected to sugar via anomeric diphosphate linkage"