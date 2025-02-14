"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation of the carboxy group of acetic acid 
    with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible acetate ester pattern
    # Acetate typically can be represented as [CX3](=O)[OX2H1] or [CX3](=O)[OX2H0] (both forms if attached)
    acetate_pattern = Chem.MolFromSmarts("C(=O)O[C,c]")

    # Define aryl structure with an -O relationship (more general than just phenol directly)
    aryl_with_o_pattern = Chem.MolFromSmarts("c[OX2H0,CX3](=O)")

    # Check if the molecule has a suitable acetate ester attached to an aryl with O
    if not mol.HasSubstructMatch(acetate_pattern):
        return False, "No acetate ester group found"

    if not mol.HasSubstructMatch(aryl_with_o_pattern):
        return False, "No aryl structure with [O] connection found"

    # Ensure the structures are linked appropriately
    # Find matches for acetate and aryl patterns (direct bond check is better than index assumption)
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    aryl_oxy_matches = mol.GetSubstructMatches(aryl_with_o_pattern)

    # Iterate through detected structures to determine connectivity
    for acetate in acetate_matches:
        acetate_carbonyl_oxygen = acetate[2]  # Indexed based on "C(=O)O[C,c]"
        acetate_oxygen = acetate[1]  # Carbonyl (O) in acetate

        for aryl in aryl_oxy_matches:
            aryl_oxygen = aryl[1]  # Linked (O) in aryl pattern

            # Check direct connectivity between acetate and aryloxy oxygen
            if mol.GetBondBetweenAtoms(acetate_oxygen, aryl_oxygen) is not None:
                return True, "Contains an appropriate linkage of an acetate ester to an aryloxy group"

    return False, "Acetate ester and aryl with oxy groups are not properly connected"

# Example usage
print(is_phenyl_acetates("CC(=O)Oc1ccccc1"))  # Should return True