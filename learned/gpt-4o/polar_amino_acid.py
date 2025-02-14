"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid contains side chains capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized amino acid backbone pattern including isotopic versions
    aa_backbone = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # Using chirality to ensure alpha-carbon is included
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "No amino acid backbone found"

    # Define specific side chains for polar amino acids
    polar_side_chains = {
        "Ser": "[CX3](O)C",    # Serine
        "Thr": "[CX4](O)C",    # Threonine
        "Asn": "[CX3](N)C(=O)", # Asparagine
        "Gln": "[CX3](N)CCC(=O)", # Glutamine
        "Cys": "[CX3](S)C",    # Cysteine
        "Tyr": "c(O)cc",       # Tyrosine (phenolic OH on aromatic ring)
        "Asp": "CC(=O)O",      # Aspartic Acid (contains COOH side group)
        "Glu": "CCC(=O)O",     # Glutamic Acid (contains COOH side group)
        "His": "Cc1c[nH]cn1",  # Histidine (contains imidazole ring)
        "Lys": "CCCC[N+]",     # Lysine (contains amine)
        "Arg": "CNC(=N)N"      # Arginine (includes complex amine)
    }

    # Check if the side chain matches any known polar side chain pattern
    for name, pattern in polar_side_chains.items():
        side_chain = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(side_chain):
            return True, f"Polar side chain '{name}' identified"

    return False, "No polar groups found in side chain"