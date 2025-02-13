"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:36976 fatty acid
Fatty acids are aliphatic monocarboxylic acids derived from or contained in esterified form in an animal or vegetable fat, oil or wax.
Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched and even-numbered), which may be saturated or unsaturated.
By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[OX2H][CX3](=O)[#6]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing carboxylic acid group"
    
    # Check for aliphatic chain
    aliphatic_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if len(aliphatic_matches) == 0:
        return False, "No aliphatic chain found"
    
    # Count carbons in the longest aliphatic chain
    chain_lengths = [len(chain) for chain in AllChem.GetUnsatRingClosures(mol, aliphatic_pattern)]
    longest_chain = max(chain_lengths)
    
    # Check chain length (4 to 28 carbons)
    if longest_chain < 4 or longest_chain > 28:
        return False, f"Aliphatic chain length {longest_chain} is outside the expected range (4-28)"
    
    # Check for cyclic structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures"
    
    # Check for branching (only linear chains)
    branched_pattern = Chem.MolFromSmarts("[CX4H2]~[CX3]")
    if mol.HasSubstructMatch(branched_pattern):
        return False, "Contains branched aliphatic chains"
    
    return True, "Molecule meets the criteria for a fatty acid"