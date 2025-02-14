"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is an amino acid not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Basic amino acid structure check: N-C-C(=O)-O or N-C-C(=O)[O-] or N-C-C(O)=O
    # Relaxed pattern to include different forms of amino acids, including beta and gamma.
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[OX2,OX1-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
       return False, "Missing basic amino acid structure"


    # 2. Check if it is a standard proteinogenic amino acid
    # Simplified exclusion list, which checks the basic substructure of proteinogenic AA
    proteinogenic_amino_acids_smarts = [
        "[NX3;H2,H1][C@H]([CH3])[CX3](=[OX1])[OX2]",  # L-Alanine
        "[NX3;H2,H1][C@@H](CC(=O)N)[CX3](=[OX1])[OX2]",  # L-Asparagine
        "[NX3;H2,H1][C@H](CC(=O)O)[CX3](=[OX1])[OX2]",  # L-Aspartic acid
        "[NX3;H2,H1][C@H](Cc1c[nH]cn1)[CX3](=[OX1])[OX2]",  # L-Histidine
        "[NX3;H2,H1][C@@H](Cc1ccccc1)[CX3](=[OX1])[OX2]",  # L-Phenylalanine
        "[NX3;H2,H1][C@H](CS)[CX3](=[OX1])[OX2]",  # L-Cysteine
        "[NX3;H2,H1][C@@H](CC(=O)N)[CX3](=[OX1])[OX2]",  # L-Glutamine
        "[NX3;H2,H1][C@H](CCC(=O)O)[CX3](=[OX1])[OX2]",  # L-Glutamic acid
        "[NX3;H2,H1][C@H](C(C)C)[CX3](=[OX1])[OX2]",  # L-Valine
        "[NX3;H2,H1][C@H](CCCNC(=N)N)[CX3](=[OX1])[OX2]", # L-Arginine
        "[NX3;H2,H1][C@H](CO)[CX3](=[OX1])[OX2]",  # L-Serine
        "[NX3;H2,H1][C@H]([C@@H](O)C)[CX3](=[OX1])[OX2]",  # L-Threonine
        "[NX3;H2,H1][C@H](CCSC)[CX3](=[OX1])[OX2]",  # L-Methionine
        "[NX3;H2,H1][C@H](C(C)CC)[CX3](=[OX1])[OX2]",  # L-Leucine
        "[NX3;H2,H1][C@H](C(C)(C)C)[CX3](=[OX1])[OX2]",  # L-Isoleucine
        "[NX3;H2,H1][C@H](Cc1c[nH]c2ccccc12)[CX3](=[OX1])[OX2]", # L-Tryptophan
        "[NX3;H2,H1][C@H](Cc1ccc(O)cc1)[CX3](=[OX1])[OX2]",  # L-Tyrosine
        "N1[C@@H]2[C@H](C[C@@H]1C(=O)O)N[CH2]2",  # L-Proline
    ]

    for smarts in proteinogenic_amino_acids_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
           return False, "It is a standard proteinogenic amino acid."
       
    return True, "It is a non-proteinogenic amino acid."