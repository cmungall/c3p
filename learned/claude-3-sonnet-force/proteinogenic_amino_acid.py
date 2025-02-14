"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:33709 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# SMARTS patterns for proteinogenic amino acids
amino_acid_patterns = {
    'glycine': '[N@H][C@@H](C(=O)O)=O',
    'alanine': '[N@H][C@@H](C)C(=O)O',
    'valine': '[N@H][C@@H](C(C)C)C(=O)O',
    'leucine': '[N@H][C@@H](CC(C)C)C(=O)O',
    'isoleucine': '[N@H][C@@H](C(C)CC)C(=O)O',
    'serine': '[N@H][C@@H](CO)C(=O)O',
    'threonine': '[N@H][C@@H](C(O)C)C(=O)O',
    'cysteine': '[N@H][C@@H](CS)C(=O)O',
    'proline': '[N@H]1C[C@@H](CC1)C(=O)O',
    'aspartic acid': '[N@H][C@@H](CC(=O)O)C(=O)O',
    'asparagine': '[N@H][C@@H](CC(N)=O)C(=O)O',
    'glutamic acid': '[N@H][C@@H](CCC(=O)O)C(=O)O',
    'glutamine': '[N@H][C@@H](CCC(N)=O)C(=O)O',
    'arginine': '[N@H][C@@H](CCCNC(N)=N)C(=O)O',
    'histidine': '[N@H][C@@H](Cc1cnc[nH]1)C(=O)O',
    'lysine': '[N@H][C@@H](CCCCN)C(=O)O',
    'phenylalanine': '[N@H][C@@H](Cc1ccccc1)C(=O)O',
    'tyrosine': '[N@H][C@@H](Cc1ccc(O)cc1)C(=O)O',
    'tryptophan': '[N@H][C@@H](Cc1c[nH]c2c1cccc2)C(=O)O',
    'methionine': '[N@H][C@@H](CCSC)C(=O)O',
    'selenocysteine': '[N@H][C@@H](C[Se]C)C(=O)O',
    'pyrrolysine': '[N@H][C@@H](CCCCNC(=O)[C@H]1CCCN1)C(=O)O',
    'N-formylmethionine': 'C[N@H]C(=O)[C@@H](CCSC)C(=O)O'
}

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule matches any of the proteinogenic amino acid patterns
    for amino_acid, pattern in amino_acid_patterns.items():
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return True, f"Molecule matches the SMARTS pattern for {amino_acid}"

    return False, "Molecule does not match any of the proteinogenic amino acid patterns"