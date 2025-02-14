"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:36370 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40) derived from psi,psi-carotene by various modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C40 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 40:
        return False, f"Molecule does not contain 40 carbons, got {c_count}"

    # Check for linear carbon backbone
    carbon_backbone = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_backbone):
        return False, "Molecule does not have a linear carbon backbone"

    # Look for cyclizations, oxidations, or rearrangements from psi,psi-carotene
    cyclized = mol.HasSubstructMatch(Chem.MolFromSmarts("[C&r5,r6]"))
    oxidized = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX1]"))
    rearranged = not mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]"))

    if not (cyclized or oxidized or rearranged):
        return False, "No modifications from psi,psi-carotene found"

    # Check for common functional groups
    has_hydroxy = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))
    has_epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts("[O;r3]"))
    has_ketone = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[!O]"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2]"))
    has_glycoside = mol.HasSubstructMatch(Chem.MolFromSmarts("C[C@@H]1[C@H]([C@@H]([C@H](O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)O)O)O)O)O)O1)O"))

    # Classify based on common carotenoid patterns
    if has_hydroxy or has_epoxide or has_ketone or has_ester or has_glycoside:
        return True, "Carotenoid skeleton with characteristic functional groups"
    else:
        return True, "Carotenoid skeleton, likely a simple carotene"