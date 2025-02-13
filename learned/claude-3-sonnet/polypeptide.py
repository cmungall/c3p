"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
A peptide containing ten or more amino acid residues.
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count amino acid residues
    amino_acids = Chem.MolFromSmarts("[N;X3](C(=O))[C@H](N)[C@H]([C,N,O])[C,N,O]")
    matches = mol.GetSubstructMatches(amino_acids)
    num_residues = len(matches)
    
    # Check if polypeptide (>= 10 amino acid residues)
    if num_residues >= 10:
        return True, f"Contains {num_residues} amino acid residues"
    else:
        return False, f"Only contains {num_residues} amino acid residues, need at least 10"

# Some key parts:
# 1. Parse the SMILES string into an RDKit molecule object
# 2. Define a SMARTS pattern to match amino acid residues
# 3. Use GetSubstructMatches to find all matches of the pattern
# 4. Count the number of matches (amino acid residues)
# 5. Classify as polypeptide if there are 10 or more residues