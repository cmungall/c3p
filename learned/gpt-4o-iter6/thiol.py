"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is defined as having a thiol group (-SH) attached to a carbon atom
    in aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): True if molecule is a thiol, False otherwise with reason
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for thiol group (-SH) directly attached to carbon
    thiol_pattern = Chem.MolFromSmarts("[CH2,CH,CH3,cH]-[S;H1]")  # Considering sp3 hybridized carbons
    
    # Check for the presence of thiol group strictly following the new pattern
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains thiol group (-SH) attached to carbon"
    
    return False, "No thiol group attached to carbon found"

# Test the function with some SMILES examples
examples = [
    ("SC(CCCCC)C", True),  # 2-Heptanethiol
    ("CO[C@@H]1\\C=C\\C=C(C)\\Cc2cc(OC)c(Cl)c(c2)N(C)C(=O)C[C@H](OC(=O)[C@H](C)N(C)C(=O)CCS)[C@]2(C)O[C@H]2[C@H](C)[C@@H]2C[C@@]1(O)NC(=O)O2", True),  # mertansine
    ("SC1=CC=C(F)C=C1", True)  # 4-Fluorothiophenol
]

for smiles, expected in examples:
    result, reason = is_thiol(smiles)
    print(f"SMILES: {smiles} --> Is Thiol: {result}, Expected: {expected}, Reason: {reason}")