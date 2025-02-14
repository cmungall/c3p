"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    Saturated fatty acids have a terminal carboxylic acid group (-COOH) and only carbon-carbon single bonds with long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the terminal carboxylic acid group. The carbon of the carboxyl group needs to be a CX3 carbon atom
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for absence of C=C or C#C bonds
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]") # checking for double bonds
    if mol.HasSubstructMatch(unsaturated_pattern):
          return False, "Contains carbon-carbon double bonds"
    
    triple_pattern = Chem.MolFromSmarts("[CX2]#[CX2]") # checking for triple bonds
    if mol.HasSubstructMatch(triple_pattern):
          return False, "Contains carbon-carbon triple bonds"
    
    # Check for at least 4 carbons in the chain using rdkit tools.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be a fatty acid"

    # Check that the structure has only one carboxylic acid group
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, must be exactly 1"
    
    # Check if there are other heteroatoms besides O and check for other functional groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 8: #Check for atoms that are not C,H or O
                return False, "Contains heteroatoms other than C, H or O"
        if atom.GetAtomicNum() == 8:  #if it's an O atom, it needs to be part of the COOH group
            is_acid_oxygen = False
            for match in acid_matches: # Check if the O atoms belong to the COOH
                acid_match_atoms = mol.GetSubstructMatch(acid_pattern) #Returns the atom indexes of the match
                if atom.GetIdx() in acid_match_atoms:
                     is_acid_oxygen = True
                     break
            if not is_acid_oxygen:
                return False, "Contains O atoms outside of carboxyl group"

    return True, "Contains a carboxylic acid group and only carbon-carbon single bonds"