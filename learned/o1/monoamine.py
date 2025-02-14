"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the monoamine SMARTS pattern
    # This matches a primary or secondary amine (NH2 or NH) not bonded to a carbonyl group,
    # connected via two aliphatic carbons (CH2-CH2) to an aromatic ring.
    monoamine_smarts = '[N;H2,H1;!$(N-C=O)]-[CH2]-[CH2]-[a]'

    # Create the SMARTS pattern molecule
    pattern = Chem.MolFromSmarts(monoamine_smarts)

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains amino group connected to aromatic ring by a two-carbon chain"
    else:
        return False, "Does not contain monoamine functional group"