"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid has a carboxylic acid group, a hydrocarbon chain,
    and only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ring structures - not allowed
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[R]")):
        return False, "Molecule contains rings"

    #Check for other atoms except C,H,O
    other_atoms = Chem.MolFromSmarts("[!#1!#6!#8]")
    if mol.HasSubstructMatch(other_atoms):
         return False, "Molecule contains non C, H, or O atoms"

    # Check for minimum chain length of at least 4 carbons excluding the carbonyl carbon (C=O).
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)
    if num_carbons < 4:
       return False, "Too few carbons to be a fatty acid"

    #Find carbon connected to carboxyl group
    chain_carbon_pattern = Chem.MolFromSmarts("CC(=O)O")
    chain_carbon_matches = mol.GetSubstructMatches(chain_carbon_pattern)
    if not chain_carbon_matches:
        return False, "Carboxylic acid not attached to carbon chain"

    # Ensure all branches are methyl branches.
    # We look for carbon atoms with more than 2 neighbors, and make sure they are connected only to CH3s, H or C.
    branch_pattern = Chem.MolFromSmarts("[CX4]([CH3])([CH3,H,C])([CH3,H,C])")
    branching_carbons_with_2methyls = mol.GetSubstructMatches(branch_pattern) # these are carbons with 2 methyls.
    branch_pattern_1methyl = Chem.MolFromSmarts("[CX4]([CH3])([CH,H])([CH,H])") # these are carbons with 1 methyl
    branching_carbons_with_1methyl = mol.GetSubstructMatches(branch_pattern_1methyl)
    branch_pattern_nomet = Chem.MolFromSmarts("[CX4]([CH,H])([CH,H])([CH,H])") #these are no methyl branch but it can be in the carbon chain.
    branching_carbons_nomet = mol.GetSubstructMatches(branch_pattern_nomet)

    #Check for non-methyl branches.
    non_methyl_branch_pattern = Chem.MolFromSmarts("[CX4]([!#1,#6])([!#1,#6,#8])([!#1,#6,#8])")
    if mol.HasSubstructMatch(non_methyl_branch_pattern):
        return False, "Molecule contains non-methyl branches"

    return True, "Methyl-branched fatty acid"