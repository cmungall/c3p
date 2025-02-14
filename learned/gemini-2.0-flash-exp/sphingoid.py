"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid has a long hydrocarbon chain, a primary amino group, and a vicinal
    hydroxyl/carbonyl group with the core structure of a C-C(OH/carbonyl)-C(N) pattern
    and sometimes with an extra OH nearby.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10 or carbon_count > 50:
       return False, f"Carbon chain length ({carbon_count}) is outside the typical range for sphingoids (10-50)"

   # Check for the presence of primary amine group
    amino_pattern = Chem.MolFromSmarts("[NX3H2]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No primary amino group found"
    
    # Check for the core sphingoid pattern: C-C(OH/carbonyl)-C(N) or related patterns
    # Defining more specific patterns with the primary amine and neighboring atoms
    # Pattern 1: C-C(OH)-C(N)-C
    sphingoid_pattern1 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3H2])[CX4]")
    # Pattern 2: C-C(=O)-C(N)-C
    sphingoid_pattern2 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3H2])[CX4]")
    # Pattern 3: C=C-C(OH)-C(N)-C
    sphingoid_pattern3 = Chem.MolFromSmarts("[CX4]=[CX3][CX4]([OH1])[CX4]([NX3H2])[CX4]")
    # Pattern 4: C=C-C(=O)-C(N)-C
    sphingoid_pattern4 = Chem.MolFromSmarts("[CX4]=[CX3][CX4](=[OX1])[CX4]([NX3H2])[CX4]")
    # Pattern 5: C-C(OH)-C(N)-C-OH
    sphingoid_pattern5 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3H2])[CX4][OH1]")
    # Pattern 6: C-C(OH)-C(N)-C(=O)
    sphingoid_pattern6 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3H2])[CX3](=[OX1])")
     # Pattern 7: C-C(OH)-C(N)-C-C-OH
    sphingoid_pattern7 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3H2])[CX4][CX4][OH1]")
    # Pattern 8: C-C-C(OH)-C(N)-C
    sphingoid_pattern8 = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH1])[CX4]([NX3H2])[CX4]")
    # Pattern 9: C-C(OH)-C(N)-C=C
    sphingoid_pattern9 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3H2])[CX3]=[CX4]")
    # Pattern 10: C-C(=O)-C(N)-C=C
    sphingoid_pattern10 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3H2])[CX3]=[CX4]")
    # Pattern 11: C-C-C(OH)-C(N)-C=C
    sphingoid_pattern11 = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH1])[CX4]([NX3H2])[CX3]=[CX4]")


    if not (mol.HasSubstructMatch(sphingoid_pattern1) or mol.HasSubstructMatch(sphingoid_pattern2) or 
            mol.HasSubstructMatch(sphingoid_pattern3) or mol.HasSubstructMatch(sphingoid_pattern4) or 
            mol.HasSubstructMatch(sphingoid_pattern5) or mol.HasSubstructMatch(sphingoid_pattern6) or 
            mol.HasSubstructMatch(sphingoid_pattern7) or mol.HasSubstructMatch(sphingoid_pattern8) or
            mol.HasSubstructMatch(sphingoid_pattern9) or mol.HasSubstructMatch(sphingoid_pattern10) or
            mol.HasSubstructMatch(sphingoid_pattern11)):
          return False, "Did not find the necessary hydroxyl or carbonyl groups near the primary amino group"
    
    # Check for a long chain with at least 10 carbons, and also some rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
         return False, "Chain too short to be a sphingoid"
    
    return True, "Contains long hydrocarbon chain, an amino group, and vicinal hydroxyl/carbonyl groups"