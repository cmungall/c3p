"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:18051 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is an amide-linked fatty acid to a sphingoid base (long-chain amino alcohol).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid base backbone pattern
    sphingoid_pattern = Chem.MolFromSmarts("[NH2][CX4;H2]([CX4])[CX4]([CX3])[CX3]=[CX3]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base backbone found"

    # Look for amide group (-C(=O)N-) linking fatty acid
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # Check for long hydrocarbon chains (fatty acid and sphingoid base)
    hydrocarbon_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    hydrocarbon_matches = mol.GetSubstructMatches(hydrocarbon_pattern)
    if len(hydrocarbon_matches) < 2:
        return False, f"Missing long hydrocarbon chains, got {len(hydrocarbon_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be ceramide"

    # Check molecular weight - ceramides typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for ceramide"

    # Count carbons and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for ceramide"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (sphingoid base)"

    return True, "Contains fatty acid amide-linked to a sphingoid base"