"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_4_hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Match flavanone structure: C6-C3-C6 with a ketone (4-carbon)
    flavanone_core_smarts = "C1=CC(=CC=C1)[C@H]2CC(=O)c3ccccc3O2"  # Core structure with stereochemistry variability
    flavanone_core_mol = Chem.MolFromSmarts(flayanone_core_smarts)
    if not mol.HasSubstructMatch(flayanone_core_mol):
        return False, "Missing core flavanone structure"

    # Look for the hydroxy group at the 4' position on B ring
    hydroxy_b_ring_smarts = "[OH]c1ccc(cc1)"  # General hydroxy position on phenyl ring
    b_ring_hydroxy_mol = Chem.MolFromSmarts(hydroxy_b_ring_smarts)
    hydroxy_matches = mol.GetSubstructMatches(b_ring_hydroxy_mol)

    # Ensure a hydroxy is specifically at the 4' position
    found_hydroxy_4prime = any(Chem.MolToSmiles(match) in list_4_hydroxy_cases for match in hydroxy_matches)
    if not found_hydroxy_4prime:
        return False, "No 4'-hydroxy group found on flavanone B ring"

    return True, "Molecule is a 4'-hydroxyflavanone"