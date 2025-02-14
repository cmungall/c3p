"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27022 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is an ubiquinone based on its SMILES string.
    An ubiquinone is a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain typically attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzoquinone core derived from 2,3-dimethoxy-5-methylbenzoquinone
    core_pattern = Chem.MolFromSmarts("C1=C(C(=O)C(=C(C1=O)[O,OX1])[O,OX1])C")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Missing benzoquinone core derived from 2,3-dimethoxy-5-methylbenzoquinone"

    # Look for polyprenoid side chain
    side_chain_pattern = Chem.MolFromSmarts("CC=C(C)CCC")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern, maxMatches=1)
    if not side_chain_matches:
        return False, "No polyprenoid side chain found"

    # Check if side chain is attached to the core
    for core_match in core_matches:
        for atom_idx in core_match:
            for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
                if bond.GetBeginAtomIdx() in side_chain_matches[0] or bond.GetEndAtomIdx() in side_chain_matches[0]:
                    return True, "Contains benzoquinone core derived from 2,3-dimethoxy-5-methylbenzoquinone with polyprenoid side chain attached"

    return False, "Polyprenoid side chain not attached to the benzoquinone core"