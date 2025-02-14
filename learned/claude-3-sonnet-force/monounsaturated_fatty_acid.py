"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: CHEBI:36195 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has one double or triple bond in the fatty acid chain and
    singly bonded carbon atoms in the rest of the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for single double or triple bond in chain
    bonds = mol.GetBonds()
    double_bonds = [b for b in bonds if b.GetBondType() == Chem.BondType.DOUBLE]
    triple_bonds = [b for b in bonds if b.GetBondType() == Chem.BondType.TRIPLE]
    if len(double_bonds) + len(triple_bonds) != 1:
        return False, "Found multiple or no double/triple bonds"

    # Check for singly bonded carbons in chain
    singly_bonded_carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalDegree() == 4]
    if len(singly_bonded_carbons) < 2:
        return False, "Insufficient singly bonded carbons in chain"

    # Check chain length (>= 4 carbons)
    chain_length = max([len(Chem.GetMolFrags(mol)[0]), len(Chem.GetMolFrags(mol)[1])])
    if chain_length < 4:
        return False, "Chain length too short for fatty acid"

    return True, "Contains one double/triple bond and singly bonded carbons in the fatty acid chain"