"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:36973 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative where one or more oxygen or hydroxy groups
    are replaced by sulfur or -SR, where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbohydrate backbone
    carbohydrate_pattern = Chem.MolFromSmarts("[C;X4;R][C;X4;R][C;X4;R]")  # Chain of 3 carbons with 4 substituents
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate backbone found"

    # Look for sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if sulfur atom is directly bonded to carbohydrate backbone
    thiosugar_pattern = Chem.MolFromSmarts("[C;X4;R][S;X2]")  # Sulfur atom attached to carbohydrate carbon
    for sulfur in sulfur_atoms:
        if mol.HasSubstructMatch(thiosugar_pattern, atomIdxList=[sulfur.GetIdx()]):
            return True, "Molecule contains a sulfur atom replacing an oxygen or hydroxy group in the carbohydrate backbone"

    return False, "Sulfur atom(s) not directly attached to the carbohydrate backbone"