"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid derivative with four double bonds,
    three of which are conjugated, derived from arachidonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a leukotriene, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a 20-carbon chain (backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Incorrect number of carbons: found {c_count}; not a leukotriene"

    # 2. Check for four double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 4:
        return False, f"Less than 4 double bonds: found {len(double_bond_matches)}; not a leukotriene"

    # 3. Check for a conjugated triene system. Allow for dihydro forms as well.
    # Look for a system of 3 conjugated double bonds plus one isolated double bond.
    # The conjugated system can have an extra single bond in the middle.
    conjugated_triene_patterns = [
        Chem.MolFromSmarts("C=C-C=C-C=C"), # strict triene
        Chem.MolFromSmarts("C=C-C=C-CC=C"), # one single bond break
        Chem.MolFromSmarts("C=CC-C=C-C=C"), # one single bond break
        Chem.MolFromSmarts("C=CC=CC=C")  # 3 conjugated bonds connected
    ]
    found_triene = False
    for pattern in conjugated_triene_patterns:
      if mol.HasSubstructMatch(pattern):
        found_triene = True
        break

    if not found_triene:
        return False, "No conjugated triene system found."


    # 4. Check for a long chain (fatty acid like)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long chain present"
        
    # 5. Check for at least one hydroxyl and one carboxylic acid
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

    if not mol.HasSubstructMatch(hydroxyl_pattern):
         return False, "No hydroxyl group found"
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Meets the criteria for a leukotriene: 20 carbons, four double bonds with 3 conjugated, a hydroxyl and a carboxylic group."