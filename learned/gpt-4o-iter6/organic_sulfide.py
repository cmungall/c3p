"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide (thioether) has the structure R-S-R', where R and R' are carbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an organic sulfide (thioether): sulfur bonded to two carbon atoms
    sulfide_pattern = Chem.MolFromSmarts("[C]-[S]-[C]")

    # Avoid classification for sulfur in sulfoxides and sulfones, which have oxygen atoms attached
    sulfoxide_sulfone_or_sulfate_pattern = Chem.MolFromSmarts("[S](=O)")

    if mol.HasSubstructMatch(sulfide_pattern):
        if mol.HasSubstructMatch(sulfoxide_sulfone_or_sulfate_pattern):
            return False, "Contains sulfur in a non-sulfide form, such as sulfoxide or sulfone"
        return True, "Contains an organic sulfide (thioether) group"

    return False, "No organic sulfide (thioether) group found"

# Example usages (selected examples illustrating different cases)
examples = [
    "CC(C)(CCO)SC[C@H](N)C(O)=O",  # felinine
    "CC(C)(O)c1ccccc1CC[C@@H](SCC1(CC1)CC([O-])=O)c1cccc(\\C=C\\c2ccc3ccc(Cl)cc3n2)c1",  # montelukast(1-)
    "CCN1C=NC2=C1N=C(N=C2NCCC3=CN=CC=C3)C#N",  # 9-ethyl-6-[2-(3-pyridinyl)ethylamino]-2-purinecarbonitrile
    "CSC1=N[C@](C)(C(=O)N1Nc1ccccc1)c1ccccc1",  # fenamidone
]

for smiles in examples:
    is_sulfide, reason = is_organic_sulfide(smiles)
    print(f"SMILES: {smiles}, is_organic_sulfide: {is_sulfide}, reason: {reason}")