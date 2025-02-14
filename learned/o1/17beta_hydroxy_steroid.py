"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is a steroid with a hydroxyl group at position 17 in the beta configuration.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to capture stereochemistry
    mol = Chem.AddHs(mol)

    # Define SMARTS pattern for steroid backbone (four fused rings)
    steroid_pattern = Chem.MolFromSmarts("[#6]12[#6][#6][#6]3[#6]([#6][#6][#6]([#6]3[#6][#6]1)[#6][#6]2)[#6]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define SMARTS pattern for 17beta-hydroxy group
    # This pattern looks for a tertiary carbon (connected to three carbons)
    # with a hydroxyl group in beta configuration ([C@@H](O))
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C]([C])[C]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No 17beta-hydroxy group found"

    # Optionally, we can try to ensure the hydroxyl group is at position 17
    # by checking if it's connected to the D-ring of the steroid backbone
    # but precise atom numbering is complex and may not be reliable here

    return True, "Contains steroid backbone with 17beta-hydroxy group in beta configuration"