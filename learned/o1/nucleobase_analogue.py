"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:25438 nucleobase analogue
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define specific SMARTS patterns for nucleobases and their analogues
    nucleobase_patterns = [
        # Purines
        Chem.MolFromSmarts("c1ncnc2ncnn12"),  # Adenine
        Chem.MolFromSmarts("c1nc2c(n1)nc(nc2=O)N"),  # Guanine
        Chem.MolFromSmarts("c1nc2c(n1)ncnc2=O"),     # Hypoxanthine
        Chem.MolFromSmarts("c1nc2c(n1)ncnc2=O"),     # Xanthine
        # Pyrimidines
        Chem.MolFromSmarts("c1cc(nc(=O)[nH]1)N"),     # Cytosine
        Chem.MolFromSmarts("c1cc(nc(=O)[nH]1)C"),     # Thymine
        Chem.MolFromSmarts("c1cc(nc(=O)[nH]1)O"),     # Uracil
        # Modified nucleobases
        Chem.MolFromSmarts("c1[nH]c(=O)c[nH]c1=O"),   # 5-fluorouracil and analogues
        Chem.MolFromSmarts("c1cc(nc(=O)[nH]1)[C,N,O]"),  # Substituted pyrimidine
        Chem.MolFromSmarts("c1ncnc2ncnc(=O)c12"),     # Purine analogues with oxygen
        Chem.MolFromSmarts("c1ncnc2ncnn12"),          # Purine analogues
    ]
    
    # Check for matches with any of the nucleobase patterns
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches nucleobase analogue pattern"
    
    return False, "Does not match nucleobase analogue patterns"

# The above code refines the SMARTS patterns to specifically match nucleobases and their analogues, avoiding overgeneralization.

# Short Summary:
# The previous attempt failed because it overgeneralized by classifying any molecule with an aromatic nitrogen-containing ring as a nucleobase analogue, leading to many false positives.
# To improve, we defined more specific SMARTS patterns that accurately represent nucleobases and their analogues.
# This approach reduces false positives by only matching molecules that have core structures similar to known nucleobases.

# The function 'is_nucleobase_analogue' checks if the input molecule matches any of the specific nucleobase patterns. If a match is found, it returns True with a reason; otherwise, it returns False.

# Finally, we avoid using overly broad criteria and focus on the structural features that are characteristic of nucleobase analogues.

# Note: This function assumes that the input SMILES string is valid and that RDKit is properly installed.

# Try running the function with some example SMILES strings:
# print(is_nucleobase_analogue("N1C=NC2=C1N=CN=C2"))  # Adenine
# print(is_nucleobase_analogue("O=C1NC(=O)NC=C1"))    # Uracil
# print(is_nucleobase_analogue("CC(=O)O"))            # Acetic acid (should return False)

# Answer Ends Here