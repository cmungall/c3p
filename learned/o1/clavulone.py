"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is an esterified prostanoid obtained from marine corals, generally characterized by a cyclopentenone ring, ester groups, long aliphatic chains, and often halogen substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentenone ring (five-membered ring with ketone)
    # Allow substitutions on ring carbons
    cyclopentenone_pattern = Chem.MolFromSmarts('O=C1CCCC1')  # Five-membered ring with ketone
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone ring found"

    # Check for ester groups (-C(=O)O- or -OC(=O)-)
    ester_pattern = Chem.MolFromSmarts('[$([CX3](=O)[OX2H0]),$([OX2H0][CX3]=O)]')  # Ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester groups found"

    # Check for long aliphatic chains (chain of at least 6 carbons)
    chain_pattern = Chem.MolFromSmarts('C' + ('C' * 5) + '[C,C+]')  # Chain of at least 6 carbons
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long aliphatic chain found"

    # Check for halogen substituents (F, Cl, Br, I) anywhere in the molecule
    halogen_pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
    if not mol.HasSubstructMatch(halogen_pattern):
        return False, "No halogen substituents found"

    # Check if molecule has at least one ring (prostanoid core)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1:
        return False, "Too few rings for a clavulone"

    return True, "Molecule matches features characteristic of clavulones"

__metadata__ = {   'chemical_class': {   'name': 'clavulone',
                              'definition': 'A class of esterified prostanoids obtained from marine corals.'}
}