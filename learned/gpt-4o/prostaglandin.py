"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are C20 compounds derived from prostanoic acid, generally featuring a
    cyclopentane ring and various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Prostaglandins usually have 20 carbons; too few carbons found"
    
    # Check for cyclopentane ring (C5 ring)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"
    
    # Check for characteristic functional groups: hydroxyl or carbonyl groups on the ring, or nearby
    oxy_groups_pattern = Chem.MolFromSmarts("[OH1,CX3](=O)")
    oxy_matches = mol.GetSubstructMatches(oxy_groups_pattern)
    if len(oxy_matches) < 1:
        return False, "Missing hydroxyl or carbonyl groups usually present in prostaglandins"
    
    # Count double bonds
    db_count = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if db_count < 2:
        return False, "Prostaglandins typically have double bonds; too few found"

    return True, "Molecule contains a C20 backbone with a cyclopentane ring and key functional groups typical of prostaglandins"