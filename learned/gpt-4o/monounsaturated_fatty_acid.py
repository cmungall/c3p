"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the carbon chain with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Find all carbon chains
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Store all bonds
    double_triple_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() in [2.0, 3.0]:
            atom1 = bond.GetBeginAtom().GetAtomicNum()
            atom2 = bond.GetEndAtom().GetAtomicNum()
            if atom1 == 6 and atom2 == 6:  # Only count C=C or C#C
                double_triple_bonds.append(bond)

    if len(double_triple_bonds) != 1:
        return False, f"Found {len(double_triple_bonds)} double/triple bonds, require exactly one"
    
    return True, "Molecule is a monounsaturated fatty acid with one double/triple bond in the main carbon chain"