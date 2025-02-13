"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:36713 saturated fatty acid
Any fatty acid containing no carbon to carbon multiple bonds. Known to produce adverse biological effects when ingested to excess.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for aliphatic carbon chain (no rings, no multiple bonds)
    aliphatic_pattern = Chem.MolFromSmarts("[C;D3]~[C;D3]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if not aliphatic_matches:
        return False, "No aliphatic carbon chain found"

    # Check for absence of multiple bonds
    if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds()):
        return False, "Contains carbon-carbon double bonds"
    if any(bond.GetBondType() == Chem.BondType.TRIPLE for bond in mol.GetBonds()):
        return False, "Contains carbon-carbon triple bonds"

    # Check carbon chain length (at least 4 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains a carboxyl group and a saturated aliphatic carbon chain"