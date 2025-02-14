"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:37412 polyprenol phosphate

A prenol phosphate resulting from the formal condensation of the terminal allylic 
hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group (P(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("P(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for prenol chain (long chain of C=C-C units)
    prenol_pattern = Chem.MolFromSmarts("[CX3]=C[CX3]C=[CX3]")
    prenol_matches = mol.GetSubstructMatches(prenol_pattern)
    if len(prenol_matches) < 3:
        return False, "Prenol chain too short or missing"

    # Check for hydroxyl group attached to prenol chain
    hydroxy_pattern = Chem.MolFromSmarts("[OX2]C[CX3]=C[CX3]C=[CX3]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No hydroxyl group attached to prenol chain"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Prenol chain too short"

    return True, "Contains a prenol chain with a terminal hydroxyl group condensed with phosphoric acid"