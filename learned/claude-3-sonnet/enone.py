"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:37858 enone
An alpha,beta-unsaturated ketone of general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= H) 
in which the C=O function is conjugated to a C=C double bond at the alpha,beta position.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha,beta-unsaturated ketone substructure
    enone_pattern = Chem.MolFromSmarts("[CD3]=[CD2][C@]=O")
    matches = mol.GetSubstructMatches(enone_pattern)
    
    if not matches:
        return False, "No alpha,beta-unsaturated ketone substructure found"
    
    # Check for conjugation
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)
        
        if c1.GetHybridization() == Chem.HybridizationType.SP2 and \
           c2.GetHybridization() == Chem.HybridizationType.SP2 and \
           c3.GetHybridization() == Chem.HybridizationType.SP2:
            return True, "Contains an alpha,beta-unsaturated ketone with conjugated C=C and C=O bonds"
    
    return False, "No conjugated alpha,beta-unsaturated ketone found"