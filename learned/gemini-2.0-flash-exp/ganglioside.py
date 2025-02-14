"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid with one or more sialic acids attached to the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Ceramide Identification
    #   - Sphingosine base (long carbon chain with a hydroxyl and an amine group with a double bond)
    #   - Fatty acid attached via amide bond
    ceramide_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]~[CX4]~[CHX4](-[OX2])~[CHX3]~[NX3]~[CX3](=[OX1])-[CX4](-[CX4])")  # Simplified ceramide pattern
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
      return False, "No ceramide structure found."

    # 2. Sialic Acid (Neu5Ac) Detection
    #   - 9-carbon skeleton with carboxyl, N-acetyl, and hydroxyl groups.
    neu5ac_pattern = Chem.MolFromSmarts("C[C](=[O])N[C]1[C]([CH]([CH]([CH]([C]([CH]([CH]([CH]1O)O)O)O)C(=O)O)O)O")
    neu5ac_matches = mol.GetSubstructMatches(neu5ac_pattern)

    if not neu5ac_matches:
        return False, "No sialic acid (Neu5Ac) found."

    #3. Check that there is a chain between the ceramide and sialic acids
    sugar_pattern = Chem.MolFromSmarts("[CX4]([OX2])-[CX4]([OX2])-[CX4]([OX2])-[OX2]-[CX4]")  # Simplified sugar pattern
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No oligosaccharide chain between ceramide and sialic acid was found."


    # Additional checks can be added
    # For example check for a minimum number of carbons and oxygens to avoid false positives.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if c_count < 30:
        return False, "Too few carbons for ganglioside."
    if o_count < 10:
        return False, "Too few oxygens for ganglioside."
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside."
    
    return True, "Contains a ceramide, an oligosaccharide, and at least one sialic acid (Neu5Ac) unit."