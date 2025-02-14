"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:36631 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a carboxylic acid with a linear or branched carbon chain,
    containing no carbon-carbon multiple bonds.

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

    # Check for carboxylic acid group
    if not any(atom.GetAtomicNum() == 8 and sum(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()) == 1 for atom in mol.GetAtoms()):
        return False, "No carboxylic acid group found"

    # Check for saturated carbon chain
    unsaturated_bonds = mol.GetBonds() if any(bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) for bond in mol.GetBonds()) else None
    if unsaturated_bonds:
        return False, "Unsaturated bonds found"

    # Check for linear or branched carbon chain
    sssr = Chem.GetSymmSSSR(mol)
    if len(sssr) > 0:
        return False, "Cyclic structures found"

    # Check for minimum chain length
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 4:
        return False, "Carbon chain too short"

    return True, "Contains a linear or branched saturated carbon chain with a carboxylic acid group"