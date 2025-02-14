"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid contains an aromatic ring directly or indirectly connected to the alpha carbon of an amino acid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring directly bonded to alpha carbon of amino acid
    # Pattern explanation:
    # [c] - aromatic carbon
    # ~ - single bond
    # [CH] - a carbon with implicit H
    # ([NX3,NX4+]) - an amine nitrogen (can be protonated or neutral)
    # C(=O)O - carboxylic acid group
    aromatic_alpha_amino_acid_pattern = Chem.MolFromSmarts("[c]~[CH]([NX3,NX4+])C(=O)O")
    if mol.HasSubstructMatch(aromatic_alpha_amino_acid_pattern):
      return True, "Contains an aromatic ring directly bonded to the alpha carbon of an amino acid moiety"

    # Check for aromatic ring connected via a linker to the alpha carbon
    aromatic_amino_acid_linker_pattern = Chem.MolFromSmarts("[c]~*-[CH]([NX3,NX4+])C(=O)O")
    if mol.HasSubstructMatch(aromatic_amino_acid_linker_pattern):
      return True, "Contains an aromatic ring connected via a linker to the alpha carbon of an amino acid moiety"

    return False, "No aromatic amino acid moiety found"