"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoserine is a glycerophosphoserine compound having an acyl substituent
    at the 1-hydroxy position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns

    # Pattern for phosphoserine group attached at sn-3 position
    phosphoserine_smarts = "[O][P](=O)([O])[O][C@@H](CN)C(=O)O"
    phosphoserine = Chem.MolFromSmarts(phosphoserine_smarts)

    # Pattern for acyl group (ester linkage) at sn-1 position
    sn1_acyl_smarts = "OCC(=O)[C,c]"
    sn1_acyl = Chem.MolFromSmarts(sn1_acyl_smarts)

    # Pattern for free hydroxyl group at sn-2 position
    sn2_oh_smarts = "OC[C@@H](O)CO[P](=O)([O])[O][C@@H](CN)C(=O)O"
    sn2_oh = Chem.MolFromSmarts(sn2_oh_smarts)

    # Check for phosphoserine group
    if not mol.HasSubstructMatch(phosphoserine):
        return False, "Phosphoserine group not found at sn-3 position"

    # Check for acyl group at sn-1 position
    if not mol.HasSubstructMatch(sn1_acyl):
        return False, "Acyl group not found at sn-1 position"

    # Check for free hydroxyl group at sn-2 position
    if not mol.HasSubstructMatch(sn2_oh):
        return False, "Free hydroxyl group not found at sn-2 position"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoserine"