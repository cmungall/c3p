"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_bile_acid_conjugate, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for steroid core (four connected rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
        
    # Check for conjugation patterns
    
    # 1. Amino acid conjugation (look for amide bond + carboxylic acid)
    amino_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4][CX3](=O)[OX2H1]")
    
    # 2. Taurine conjugation
    taurine_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4][CX4][SX4](=O)(=O)[OX2H1,OX1-]")
    
    # 3. Sulfate conjugation
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX2H1,OX1-]")
    
    # 4. Glucuronic acid conjugation (simplified pattern)
    glucuronic_pattern = Chem.MolFromSmarts("[OX2]C1[OX2][CX4][CX4][CX4]([CX3](=O)[OX2H1])[CX4]1")
    
    # 5. Sugar conjugation (simplified pattern for glucose)
    sugar_pattern = Chem.MolFromSmarts("[OX2]C1[OX2][CX4][CX4][CX4][CX4]1")
    
    conjugation_patterns = [
        (amino_acid_pattern, "amino acid"),
        (taurine_pattern, "taurine"),
        (sulfate_pattern, "sulfate"),
        (glucuronic_pattern, "glucuronic acid"),
        (sugar_pattern, "sugar")
    ]
    
    found_conjugations = []
    for pattern, name in conjugation_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_conjugations.append(name)
            
    if not found_conjugations:
        return False, "No conjugation patterns found"
        
    # Check for hydroxyl groups (common in bile acids)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 1:
        return False, "No hydroxyl groups found on steroid core"
        
    return True, f"Bile acid conjugated with {', '.join(found_conjugations)}"