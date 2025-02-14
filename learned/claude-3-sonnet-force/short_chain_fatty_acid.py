"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:35834 short-chain fatty acid
A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6.
If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for aliphatic chain (no rings, no heteroatoms except O)
    aliphatic_pattern = Chem.MolFromSmarts("[CH2,CH3,CH,C](C)(C)(C)")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    non_aliphatic_atoms = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() not in ['C', 'H', 'O']]
    if not aliphatic_matches or non_aliphatic_atoms:
        return False, "Molecule is not aliphatic"
    
    # Count carbon chain length
    carbon_chain = max([len(Chem.FragmentOnSomething(mol, atom).GetAtoms()) for atom in aliphatic_matches])
    if carbon_chain >= 6:
        return False, "Carbon chain is too long (>= C6)"
    
    # Check for non-hydrocarbon substituents
    non_h_substituents = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and sum(n.GetSymbol() != 'H' for n in atom.GetNeighbors()) > 2]
    if non_h_substituents:
        return False, "Molecule contains non-hydrocarbon substituents"
    
    return True, "Aliphatic monocarboxylic acid with chain length < C6 and no non-hydrocarbon substituents"