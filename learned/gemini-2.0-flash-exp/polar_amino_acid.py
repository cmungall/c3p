"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid has a side chain capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True if molecule is a polar amino acid, False otherwise
                         and the reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the amino acid backbone pattern
    backbone_pattern = Chem.MolFromSmarts("[NX3;H2,H1]-[CX4]-[CX3](=[OX1])-[OX2;H1]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Molecule does not have the amino acid backbone"

    # Check for polar side chain patterns
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[OX2;H1]"),  # hydroxyl (-OH)
        Chem.MolFromSmarts("[NX3;H2,H1]"),  # amino (-NH2)
        Chem.MolFromSmarts("C(=[OX1])N"),  # amide (-CONH2)
        Chem.MolFromSmarts("C(=[OX1])[OX2;H1]"),  # carboxylic acid (-COOH)
        Chem.MolFromSmarts("[SX2;H1]"), #thiol (-SH)
        Chem.MolFromSmarts("c[nH]c"),   #imidazole
        Chem.MolFromSmarts("c1ccccc1[OH]"), #phenol
        Chem.MolFromSmarts("c1ccncc1N") #pyridine
    ]
    
    has_polar_group = False
    for pattern in polar_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            has_polar_group = True
            break
    
    if not has_polar_group:
        return False, "No polar side chain capable of forming hydrogen bonds found."


    return True, "Molecule has amino acid backbone with polar side chain"