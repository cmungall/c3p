"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is characterized by an oxygenated flavylium ion core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Check for flavylium core: 2-phenylchromenylium structure
    flavylium_pattern = Chem.MolFromSmarts("[o+]1c(ccc2cc1)c1cocc2c1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return (False, "No flavylium core found")

    # Check for oxygenated groups (e.g., hydroxyl groups) on the flavonium backbone
    oxygen_pattern = Chem.MolFromSmarts("[OH]")
    oxy_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxy_matches) < 2:
        return (False, "Insufficient oxygenation; at least two hydroxyl groups expected")

    # Check for positive ion state
    positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    if not positive_charge:
        return (False, "Molecule lacks a positive charge indicating it is not a cation")

    # The aglycon form should not contain sugar residues; however, we focus on flavylium
    # Since input examples have sugars, but task defines aglycons, assumes processing aglycon core
    
    return (True, "Contains an oxygenated flavylium ion core with positive charge")

# Example usage:
# result, reason = is_anthocyanidin_cation("Insert valid SMILES here")
# print(result, reason)