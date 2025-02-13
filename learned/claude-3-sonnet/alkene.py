"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is cyclic
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings - must be acyclic"

    # Check if molecule contains only C and H atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:  # 1=H, 6=C
            return False, "Contains non-hydrocarbon atoms"

    # Count double bonds between carbons
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds == 0:
        return False, "No carbon-carbon double bonds found"
    elif n_double_bonds > 1:
        return False, "Contains multiple double bonds"

    # Count carbons and hydrogens to verify CnH2n formula
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Calculate total H count (explicit + implicit)
    total_h_count = sum(atom.GetTotalNumHs() + atom.GetNumExplicitHs() 
                       for atom in mol.GetAtoms())
    
    # Check if follows CnH2n formula
    if total_h_count != 2 * c_count:
        return False, f"Does not follow CnH2n formula: C{c_count}H{total_h_count}"

    return True, f"Acyclic hydrocarbon with one C=C double bond, follows C{c_count}H{total_h_count} formula"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:32878',
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one '
                     'carbon-carbon double bond and the general formula CnH2n. '
                     'Acyclic branched or unbranched hydrocarbons having more than '
                     'one double bond are alkadienes, alkatrienes, etc.',
        'parents': ['CHEBI:33641']
    }
}