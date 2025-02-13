"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:28051 catechins
Catechins are members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin contains a flavan-3-ol core and various hydroxyl and methoxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavan-3-ol core
    flavan_3_ol_pattern = Chem.MolFromSmarts("[CR1]1=CC=C2C(=C1)OC(C3=CC=C(O)C=C3)C[C@@H]2[C@@H](O)[C@@H](O)CC4=CC=CC=C4")
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol core found"
    
    # Check for common catechin substituents
    allow_subs = ["c", "cO", "cOC", "c(O)c", "c(OC)c", "cc(O)", "cc(OC)"]
    sub_smarts = Chem.MolFromSmarts(f"[{','.join(allow_subs)}]")
    sub_matches = mol.GetSubstructMatches(sub_smarts)
    
    # Require at least 2 hydroxy/methoxy substituents
    if len(sub_matches) < 2:
        return False, "Not enough hydroxy/methoxy substituents"
    
    # Check for common non-catechin substructures
    ban_smarts = ["c1ccccc1", "c1nccnc1", "C#N", "C=O", "N", "S", "P", "F", "Cl", "Br", "I"]
    for smarts in ban_smarts:
        ban_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(ban_pattern):
            return False, f"Contains banned substructure: {smarts}"
        
    return True, "Contains flavan-3-ol core with hydroxy/methoxy substituents"