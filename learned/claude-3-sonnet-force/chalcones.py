"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:38236 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcone(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone with the structure ArCH=CH(=O)Ar, where Ar represents an aryl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chalcone backbone pattern (two aryl groups separated by ketone and alkene)
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)c1ccccc1")
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "No chalcone backbone found"

    # Check for aromaticity and specific atom counts
    aromatic_rings = AllChem.GetSSSR(mol)
    if len(aromatic_rings) < 2:
        return False, "Not enough aromatic rings for chalcone"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 15 or o_count != 1:
        return False, "Incorrect carbon or oxygen count"

    # Count rotatable bonds - chalcones typically have 2-3
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2 or n_rotatable > 3:
        return False, "Incorrect number of rotatable bonds"

    return True, "Contains the chalcone backbone ArCH=CH(=O)Ar"