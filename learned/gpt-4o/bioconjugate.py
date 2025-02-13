"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two distinct biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a bioconjugate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define a range of known biological moiety patterns for identification
    patterns = {
        "amino_acid": Chem.MolFromSmarts("N[C@@H](C)C(=O)O"),                    # Simplified amino acid backbones
        "fatty_acid": Chem.MolFromSmarts("C(=O)[CX4]C"),                         # General fatty acid residue
        "nucleotide": Chem.MolFromSmarts("n1cnc2c1ncnc2"),                       # DNA/RNA bases
        "sugar": Chem.MolFromSmarts("C(O)C(O)C(O)C(O)"),                         # Simple sugar backbone
        "peptide_bond": Chem.MolFromSmarts("C(=O)N"),                            # Peptide linkage
        "thiol": Chem.MolFromSmarts("CSC(C)"),                                   # Cysteine related
        "phosphate": Chem.MolFromSmarts("P(=O)(O)O")                             # Phosphate group
    }

    # Track identified biological moieties
    identified_moieties = set()

    # Check for the presence of each biological moiety pattern
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            identified_moieties.add(name)

    # Consider it a bioconjugate if two or more distinct moieties are found
    if len(identified_moieties) >= 2:
        return True, f"Contains at least two different biological moieties covalently linked: {', '.join(identified_moieties)}"
    
    return False, "Does not contain at least two distinct biological components"

# Example use (replace SMILES with an actual string to test):
# result, reason = is_bioconjugate("your_smiles_here")