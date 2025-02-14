"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined by a ketal carbon that is the only common atom of two distinct rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a spiroketal: A carbon single-bonded to two oxygens each being part of distinct rings
    spiroketal_pattern = Chem.MolFromSmarts("[C]([O])[O]")  # This needs further refinement

    # Check for spiroketal substructure in the molecule
    if not spiroketal_pattern or not mol.HasSubstructMatch(spiroketal_pattern):
        return False, "No spiroketal structure found"

    # Verify that this setup forms a valid spiro configuration with two distinct rings
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    for match in matches:
        ketal_carbon = match[0]
        connected_oxygens = [match[1], match[2]]
        rings = [set(ring) for ring in Chem.GetSymmSSSR(mol) if any(o in ring for o in connected_oxygens)]
        if len(rings) != 2 or len(rings[0].intersection(rings[1])) != 1:
            continue  # This is not a spiroketal

        return True, "Contains ketal carbon common to two distinct rings (spiroketal structure)"

    return False, "No valid spiroketal: ketal carbon is not in two distinct rings"