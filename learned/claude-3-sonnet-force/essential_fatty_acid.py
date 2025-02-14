"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:33567 essential fatty acid

Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid with an absolute dietary requirement.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon chain
    longest_chain = AllChem.FindLongestChain(mol)
    if len(longest_chain) < 12:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for multiple double bonds (polyunsaturated)
    num_double_bonds = sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    if num_double_bonds < 2:
        return False, "Fewer than 2 double bonds, not polyunsaturated"
    
    # Check for cis configuration of double bonds
    smarts_cis = Chem.MolFromSmarts("/C=C/")
    cis_matches = mol.GetSubstructMatches(smarts_cis)
    if not cis_matches:
        return False, "Double bonds not in cis configuration"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 16:
        return False, "Too few carbons for essential fatty acid"
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (carboxylic acid)"
    
    # Check molecular weight - essential fatty acids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for essential fatty acid"
    
    return True, "Polyunsaturated fatty acid with multiple cis double bonds and carboxylic acid group"