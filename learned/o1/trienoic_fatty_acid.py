"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a fatty acid with a carboxylic acid group and exactly three carbon-carbon double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)OH or -C(=O)O-)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,H]")  # Matches both -COOH and -COO(-)
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon-carbon double bonds
    num_double_bonds = 0
    for bond in mol.GetBonds():
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and
            bond.GetBeginAtom().GetAtomicNum() == 6 and
            bond.GetEndAtom().GetAtomicNum() == 6):
            num_double_bonds += 1
            
    if num_double_bonds != 3:
        return False, f"Contains {num_double_bonds} carbon-carbon double bonds, needs exactly 3"

    # Optional: Check if the hydrocarbon chain is sufficiently long (typical for fatty acids)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Contains {c_count} carbon atoms, which is too short for a typical fatty acid chain"
    
    return True, "Molecule is a trienoic fatty acid with a carboxylic acid group and three C=C double bonds"