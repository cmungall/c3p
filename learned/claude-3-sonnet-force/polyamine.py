"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:35586 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic amino compound containing two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of primary, secondary, and tertiary amino groups
    primary_amino_pattern = Chem.MolFromSmarts("[NH2;!$(NC=O)]")
    secondary_amino_pattern = Chem.MolFromSmarts("[NH1;!$(NC=O)]")
    tertiary_amino_pattern = Chem.MolFromSmarts("[NH0;!$(N([C;!$(C=[O,N,P,S])])(C))]")
    amino_count = len(mol.GetSubstructMatches(primary_amino_pattern)) + \
                  len(mol.GetSubstructMatches(secondary_amino_pattern)) + \
                  len(mol.GetSubstructMatches(tertiary_amino_pattern))

    # Check for at least two amino groups
    if amino_count >= 2:
        # Exclude molecules containing amide groups
        amide_pattern = Chem.MolFromSmarts("C(=O)N")
        if mol.HasSubstructMatch(amide_pattern):
            return False, "Contains amide group(s), not a polyamine"
        else:
            return True, f"Found {amino_count} amino groups"
    else:
        return False, f"Found only {amino_count} amino group(s), need at least 2"