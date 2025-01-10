"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins have a beta-lactam structure fused with a six-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more refined beta-lactam ring SMARTS pattern
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N[C@H]1")  # Focusing on the beta-lactam core more directly
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Define a six-membered dihydrothiazine-containing pattern typical for cephalosporins.
    # The C@ atom indicates a potential chiral center present in these types.
    dihydrothiazine_pattern = Chem.MolFromSmarts("C2[S][C@H]1N2C(=O)C=C1")  # Focusing on the six-membered thiazine ring
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found"

    return True, "Contains beta-lactam and dihydrothiazine rings characteristic of cephalosporins"

# Example testing
smiles_examples = [
    "C=1(N2[C@](SCC1/C=C/3\CCN(C3=O)[C@]4(CN(CC4)C(OCC=5OC(OC5C)=O)=O)[H])([C@@](C2=O)(NC(/C(/C6=NSC(=N6)N)=N\O)=O)[H])[H])C(O)=O",
    "CC1=C(C(=NN1CC(=O)NC2C3N(C2=O)C(=C(CS3)CSC4=NN=C(S4)C)C(=O)O)C)Cl",
    "S1[C@]2(N(C(=O)[C@H]2NC(=O)\C(=N\OC)\C=3N=C(SC3)N)C(=C(C1)CO)C(O)=O)[H]"
]

for smiles in smiles_examples:
    is_cef, reason = is_cephalosporin(smiles)
    print(f"SMILES: {smiles} - Is Cephalosporin? {is_cef}. Reason: {reason}")