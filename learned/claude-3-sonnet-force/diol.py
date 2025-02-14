"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:16549 diol
A compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
Aliphatic diols are also called glycols.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))

    # Check for specific diol substructures
    diol_patterns = [
        Chem.MolFromSmarts("[OX2H][CX4][CX4][OX2H]"),  # 1,2-diol
        Chem.MolFromSmarts("[OX2H]C1CCCCC1[OX2H]"),  # cyclic diol
        Chem.MolFromSmarts("[OX2H][CX4H2][OX2H]"),  # geminal diol
        Chem.MolFromSmarts("[OX2H][CX4][CX4][CX4][OX2H]"),  # 1,3-diol
        Chem.MolFromSmarts("[OX2H][CX4][CX4][CX4][CX4][OX2H]"),  # 1,4-diol
        # Add more diol patterns as needed
    ]

    has_diol_substructure = any(mol.HasSubstructMatch(p) for p in diol_patterns)

    # Check for disqualifying substructures
    disqualifying_patterns = [
        Chem.MolFromSmarts("[C$(C(=O)O)]=O"),  # Carboxylic acid
        Chem.MolFromSmarts("C(=O)O[CX4]"),  # Ester
        # Add more disqualifying patterns as needed
    ]

    has_disqualifying_substructure = any(mol.HasSubstructMatch(p) for p in disqualifying_patterns)

    # Make the final decision
    if hydroxy_count == 2 and has_diol_substructure and not has_disqualifying_substructure:
        return True, "Contains two hydroxy groups in a typical diol arrangement"
    elif hydroxy_count > 2:
        return False, "Contains more than two hydroxy groups, not a typical diol"
    else:
        return False, "Does not contain two hydroxy groups in a typical diol arrangement"