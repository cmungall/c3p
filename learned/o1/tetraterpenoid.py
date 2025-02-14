"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26976 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is derived from a tetraterpene (C40 backbone), possibly modified
    by rearrangement or removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbons ({c_count}) for a tetraterpenoid"

    # Count number of conjugated double bonds
    # We will approximate conjugation by counting alternating single and double bonds in a chain
    # This is a simplification and may not capture all cases
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_count += 1

    if double_bond_count < 10:
        return False, f"Too few double bonds ({double_bond_count}) for a tetraterpenoid"

    # Estimate presence of extended conjugation by calculating the fraction of double bonds
    total_bonds = mol.GetNumBonds()
    if double_bond_count / total_bonds < 0.2:
        return False, "Insufficient fraction of double bonds for extended conjugation"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a tetraterpenoid"

    # Check for the presence of isoprene units (C5H8)
    # This is challenging due to possible rearrangements, so we will attempt to find repeating units
    isoprene_smarts = Chem.MolFromSmarts("C(C)=C-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_smarts)
    if len(isoprene_matches) < 4:
        return False, "Insufficient isoprene-like units found"

    return True, "Molecule matches criteria for a tetraterpenoid (high carbon count, extended conjugation, isoprenoid pattern)"