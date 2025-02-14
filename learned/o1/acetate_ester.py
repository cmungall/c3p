"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an ester group: carbonyl carbon single bonded to an oxygen
    ester_pattern = Chem.MolFromSmarts("C(=O)[O][#0,#6,#7,#8,#9,#15,#16,#17,#35,#53]")
    if ester_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for matches of the ester pattern in the molecule
    matches = mol.GetSubstructMatches(ester_pattern)
    if not matches:
        return False, "No ester groups found"

    acetate_count = 0  # Counter for acetate ester groups

    for match in matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[2]

        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        ester_o = mol.GetAtomWithIdx(ester_o_idx)

        # Get neighbors of the carbonyl carbon excluding oxygens
        neighbors = [atom for atom in carbonyl_c.GetNeighbors() if atom.GetAtomicNum() != 8]
        
        # Exclude double-bonded oxygen (carbonyl oxygen)
        neighbors = [atom for atom in neighbors if atom.GetIdx() != ester_o_idx]

        if len(neighbors) != 1:
            continue  # Not an acetate ester if more or fewer than one neighbor besides oxygens

        alpha_c = neighbors[0]

        # Check if the alpha carbon is a methyl group (degree 1 carbon)
        if alpha_c.GetAtomicNum() == 6 and alpha_c.GetDegree() == 1:
            acetate_count += 1

    if acetate_count > 0:
        return True, f"Contains {acetate_count} acetate ester group(s)"
    else:
        return False, "No acetate ester groups found"