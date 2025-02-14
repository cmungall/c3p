"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide SMARTS pattern
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2]-[CX4H3]")

    # Check for substructure match of a methyl sulfide
    matches = mol.GetSubstructMatches(methyl_sulfide_pattern)
    if not matches:
        return False, "Molecule does not contain a methyl sulfide group"

    # Check for aliphaticity around the methyl sulfide.
    for match in matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        methyl_carbon_atom = mol.GetAtomWithIdx(match[1])

        if not sulfur_atom.GetIsAromatic() and not methyl_carbon_atom.GetIsAromatic():
            return True, "Molecule contains an aliphatic methyl sulfide group"
        
    return False, "Molecule does not contain an aliphatic methyl sulfide group"