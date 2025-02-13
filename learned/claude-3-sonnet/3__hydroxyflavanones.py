"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
Compounds with a flavanone core structure and a hydroxy substituent at position 3' of the phenyl ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavanone core structure with 3'-OH
    # [OH1] ensures we match a hydroxy group
    # The position numbering ensures we match the 3' position specifically
    flavanone_3_oh_pattern = Chem.MolFromSmarts("""
        [OH1]c1cc([C,c]2CC(=O)c3ccccc3O2)ccc1
    """)
    
    # Alternative pattern for when the OH is written in reverse order
    flavanone_3_oh_pattern_alt = Chem.MolFromSmarts("""
        c1cc([C,c]2CC(=O)c3ccccc3O2)ccc1[OH1]
    """)

    # Check for matches
    if mol.HasSubstructMatch(flavanone_3_oh_pattern) or mol.HasSubstructMatch(flavanone_3_oh_pattern_alt):
        # Verify basic flavanone structure (2 rings with ketone)
        flavanone_core = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc21")
        if not mol.HasSubstructMatch(flavanone_core):
            return False, "Has 3'-OH but missing flavanone core structure"
        
        return True, "Contains flavanone core with 3'-hydroxy substituent"

    # Check if it has flavanone core but missing 3'-OH
    if mol.HasSubstructMatch(flavanone_core):
        return False, "Has flavanone core but missing 3'-hydroxy substituent"
    
    return False, "Does not contain flavanone core structure"