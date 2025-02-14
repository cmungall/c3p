"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bond pattern (O=C-N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)[NH]")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"
    
    # Check for overall structural complexity - proxy using heavy atom count
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count < 30:
        return False, f"Too few heavy atoms ({heavy_atom_count}) for typical peptide antibiotic"
    
    # Additional feature patterns can be included here as identified
    # For simplicity, return True based on basic checks, but real-world would require complex patterns

    return True, "Contains peptide bonds with a complex structure typical of peptide antibiotics"

# Note: Real detailed classification would require specific functional group and motif checks