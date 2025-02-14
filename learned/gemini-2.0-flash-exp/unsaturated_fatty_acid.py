"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a carboxylic acid.
    It also must be a long, linear chain, with no aromatic rings and no phosphorus.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

     # Check for phosphorus atoms
    has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if has_phosphorus:
        return False, "Contains phosphorus, not a fatty acid"

    # Check for aromatic rings
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings, not a fatty acid"

    # Check for carboxylic acid group (-C(=O)O or -C(=O)[O-]) at the end of a chain
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2;!$([OX2](C)(C))!H0]")
    if not mol.HasSubstructMatch(acid_pattern):
       return False, "Does not contain a carboxylic acid group"

    # Define SMARTS patterns for C=C and C#C bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if not (mol.HasSubstructMatch(double_bond_pattern) or mol.HasSubstructMatch(triple_bond_pattern)):
        return False, "Does not contain any C=C or C#C bond"

    # Check for a long carbon chain (at least 4 carbons) outside of a ring. Also ensure that it's a linear chain.
    chain_pattern = Chem.MolFromSmarts("[!R;CX4,CX3]~[!R;CX4,CX3]~[!R;CX4,CX3]~[!R;CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
         return False, "Does not contain a sufficiently long linear carbon chain"

    # Check number of carbons
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 4:
        return False, "Too few carbon atoms"

    # If all checks pass
    return True, "Contains a carboxylic acid group, a long chain and at least one C=C or C#C bond"