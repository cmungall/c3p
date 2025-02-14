"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an arylalkylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one amino group (primary, secondary, tertiary, or quarternary)
    # and explicitly exclude amino groups involved in ring systems:
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N=*)]")  # Primary, secondary, tertiary
    ammonium_pattern = Chem.MolFromSmarts("[NX4;!$(N=*)]")  # Quarternary ammonium
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    ammonium_matches = mol.GetSubstructMatches(ammonium_pattern)

    if len(amino_matches) + len(ammonium_matches) != 1:
        return False, f"Molecule must have exactly one amino group. Found {len(amino_matches) + len(ammonium_matches)}"

    # Check for the presence of at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[c]")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
       return False, "Molecule must contain at least one aromatic ring."
       
    # Check for the two-carbon chain directly linking the amino group to an aromatic ring.
    # The chain is explicitly required to be two single-bonded saturated carbons directly connected to an aromatic atom and the N
    two_carbon_chain_pattern = Chem.MolFromSmarts("[NX3,NX4][CH2][CH2][c]")
    two_carbon_matches = mol.GetSubstructMatches(two_carbon_chain_pattern)

    if len(two_carbon_matches) != 1:
        return False, f"Molecule must have exactly one amino group connected to an aromatic ring by a two-carbon chain. Found {len(two_carbon_matches)}"
    
    # Check if any nitrogen is part of a ring.
    nitrogen_in_ring = Chem.MolFromSmarts("[NX3,NX4]1[#6,#7,#8,#16]~*1")
    for match in amino_matches + ammonium_matches: # Check for each nitrogen individually.
        nitrogen_atom = mol.GetAtomWithIdx(match[0])
        if nitrogen_atom.IsInRing():
           return False, "Amino group cannot be part of a ring system."
    
    return True, "Molecule contains one amino group attached to an aromatic ring via a two-carbon chain"