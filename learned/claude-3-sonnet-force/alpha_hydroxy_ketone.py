"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35708 alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone contains a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use a SMARTS pattern to match alpha-hydroxy ketones
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])(C)[C@H](O)")
    matches = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)

    if matches:
        # Check molecular weight as an additional filter
        mol_wt = Descriptors.MolWt(mol)
        if mol_wt < 100 or mol_wt > 1000:
            return False, "Molecular weight outside typical range for alpha-hydroxy ketones"

        # Check for presence of additional functional groups
        # (This is just an example, you may want to modify or remove this check)
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3]")):
            return False, "Contains additional nitrogen-containing functional groups"

        return True, "Contains a hydroxy group on the alpha-carbon relative to the C=O group"

    return False, "Does not contain an alpha-hydroxy ketone substructure"