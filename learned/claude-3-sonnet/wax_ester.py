"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: wax ester
A fatty acid ester resulting from the condensation of the carboxy group of a fatty acid 
with the alcoholic hydroxy group of a fatty alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for elements other than C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, and O"

    # Count oxygens - should be exactly 2 for a wax ester
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (one ester group)"

    # Look for fatty acid chain pattern on carbonyl side
    # At least 7 carbons in a chain attached to the carbonyl
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain found"

    # Look for fatty alcohol chain pattern on alcohol side
    # At least 7 carbons in a chain attached to the ester oxygen
    fatty_alcohol_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(fatty_alcohol_pattern):
        return False, "No fatty alcohol chain found"

    # Count carbons - wax esters typically have at least 16 total carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Total carbon count too low for wax ester"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for wax ester"

    # Look for problematic groups that would make it not a wax ester
    problematic_groups = [
        "[OH]",  # alcohol (except the one in ester)
        "[CX3](=O)[OX2H]",  # carboxylic acid
        "[NX3]",  # amine
        "[SX2]",  # thiol/sulfide
    ]
    
    for smart in problematic_groups:
        pattern = Chem.MolFromSmarts(smart)
        if mol.HasSubstructMatch(pattern):
            return False, "Contains additional functional groups"

    # Optional: Allow for unsaturation in the chains
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_unsaturation = mol.HasSubstructMatch(double_bond_pattern)

    return True, f"Contains one ester group connecting two long hydrocarbon chains{' with unsaturation' if has_unsaturation else ''}"