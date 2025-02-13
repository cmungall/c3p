"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as having a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update SMARTS pattern for thiol groups to be more inclusive of different alkanethiol contexts
    thiol_pattern = Chem.MolFromSmarts("[S&H1]-[C]")  # SH must be attached to any carbon, not limited to -CH2-

    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol group (-SH) attached to an alkyl group found"

    # Exclude certain complex structures that are not simple alkanethiols
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            connected_atoms = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            # Check for forbidden patterns
            if 7 in connected_atoms or 8 in connected_atoms or 15 in connected_atoms:
                return False, "Sulfur atom bonded with non-alkyl group components or complex structures, possibly involving peptide, amide, or phosphorous compounds"

    return True, "Contains a sulfanyl group, -SH, attached to an alkyl chain"