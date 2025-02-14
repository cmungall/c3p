"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18357 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine contains a catechol moiety (benzene ring with adjacent hydroxyl groups)
    and an aminoethyl side chain attached to the benzene ring, including derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol moiety pattern (benzene ring with two adjacent hydroxyl groups)
    # Allow for substitutions on the ring
    catechol_pattern = Chem.MolFromSmarts('c1cc([OH])c([OH])cc1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found (benzene ring with adjacent hydroxyl groups)"
    
    # Define aminoethyl side chain pattern attached to the aromatic ring
    # Allow for substitutions on the amino group and chain
    aminoethyl_pattern = Chem.MolFromSmarts('[$(c-[CH2]-[CH2]-[N])]')  # Aromatic carbon connected to ethylamine group
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain attached to aromatic ring found"
    
    # Ensure that the catechol moiety and aminoethyl side chain are connected to the same ring
    # Get atom indices of catechol hydroxyl oxygens
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)
    
    # Check if the aromatic carbons in catechol and aminoethyl patterns overlap
    for cat_match in catechol_matches:
        # Extract the aromatic carbons in the catechol moiety
        catechol_ring = [idx for idx in cat_match if mol.GetAtomWithIdx(idx).IsInRing()]
        for amine_match in aminoethyl_matches:
            # The first atom in the aminoethyl pattern is the aromatic carbon
            aromatic_carbon_idx = amine_match[0]
            if aromatic_carbon_idx in catechol_ring:
                return True, "Contains catechol moiety with aminoethyl side chain attached"
    
    return False, "Catechol moiety and aminoethyl side chain not properly connected"