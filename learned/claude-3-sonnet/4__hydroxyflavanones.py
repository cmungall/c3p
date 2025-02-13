"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4'
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_4_hydroxyflavanone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic flavanone core structure
    # [#6]1~[#6]~[#6](=[O])~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[O]1
    flavanone_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6](=[O])~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2~[O]1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core structure found"
    
    # Check for ketone group at position 4
    ketone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group at position 4"
    
    # Check for 4'-hydroxy group
    # Match pattern for flavanone with 4'-hydroxy specifically
    hydroxyflavanone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccc(O)cc2)Oc2ccccc12")
    if not mol.HasSubstructMatch(hydroxyflavanone_pattern):
        return False, "Missing hydroxy group at 4' position"
    
    # Additional check to ensure saturated C2-C3 bond (distinguishing from flavones)
    # The pattern matches the C2-C3 single bond specifically
    c2c3_pattern = Chem.MolFromSmarts("O=C1C[CH2][CH]([c2]ccccc2)Oc2ccccc12")
    if not mol.HasSubstructMatch(c2c3_pattern):
        return False, "C2-C3 bond must be saturated (not a flavone)"
    
    # If all checks pass, it's a 4'-hydroxyflavanone
    return True, "Molecule contains flavanone core with 4'-hydroxy group"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        ("O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1", True),  # taxifolin
        ("Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1", True),  # naringenin
        ("COc1cc(O)cc2OC(CC(=O)c12)c1ccccc1", False),  # missing 4'-OH
        ("O=C1CC(c2ccc(O)cc2)Oc2ccccc12", True),  # basic 4'-hydroxyflavanone
    ]
    
    for smiles, expected in examples:
        result, reason = is_4__hydroxyflavanones(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()