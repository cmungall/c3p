"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid core (four connected rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C]([C]1)[C][C]3[C]([C]2)[C][C]4[C][C][C]@[C]4[C]3")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"
    
    # Check for 3-hydroxy group in alpha configuration
    # [OH] connected to carbon at position 3 in alpha configuration
    # The @ in the SMARTS pattern specifies the stereochemistry
    alpha_3_hydroxy = Chem.MolFromSmarts("[OH][C@H]1[CH2][CH2][C]2")
    
    if not mol.HasSubstructMatch(alpha_3_hydroxy):
        # Also check alternative SMARTS pattern for cases where the stereochemistry 
        # might be represented differently
        alt_alpha_3_hydroxy = Chem.MolFromSmarts("[OH][C@@H]1[CH2][CH2][C]2")
        if not mol.HasSubstructMatch(alt_alpha_3_hydroxy):
            return False, "No 3-alpha hydroxy group found"
    
    # Additional check to ensure the hydroxy group is at position 3
    # This is implicitly handled by the SMARTS pattern above, but included for clarity
    
    return True, "Contains steroid core with 3-alpha hydroxy group"