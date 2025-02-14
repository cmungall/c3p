"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is any fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Count number of carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 4:
        return False, "Too few carbon atoms for a fatty acid"

    # Check for C=C double bonds
    double_bond = Chem.MolFromSmarts('C=C')
    has_double_bond = mol.HasSubstructMatch(double_bond)

    # Check for C#C triple bonds
    triple_bond = Chem.MolFromSmarts('C#C')
    has_triple_bond = mol.HasSubstructMatch(triple_bond)

    if not (has_double_bond or has_triple_bond):
        return False, "No carbon-carbon double or triple bonds found"

    return True, "Molecule is an unsaturated fatty acid"