"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: CHEBI:51766 long-chain fatty alcohol
A fatty alcohol with a chain length ranging from C13 to C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with an -OH group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an -OH group
    oh_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxy (-OH) group found"

    # Count carbon atoms in the longest chain
    sssr = Chem.GetSSSR(mol)
    longest_chain_length = max(len(ring) for ring in sssr)
    if not (13 <= longest_chain_length <= 22):
        return False, f"Longest chain length is {longest_chain_length}, not in the range C13 to C22"

    # Check if longest chain is aliphatic (no ring atoms)
    longest_chain = max(sssr, key=len)
    if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in longest_chain):
        return False, "Longest chain contains ring atoms"

    # Check for other functional groups (e.g. esters, ketones)
    allowed_patterns = [Chem.MolFromSmarts("[OX1H]"), Chem.MolFromSmarts("[CX4H3]")]
    for atom in mol.GetAtoms():
        if not any(atom.HasMatch(pattern) for pattern in allowed_patterns):
            return False, "Contains unexpected functional groups"

    return True, "Molecule contains a straight-chain aliphatic alcohol with a chain length between C13 and C22"