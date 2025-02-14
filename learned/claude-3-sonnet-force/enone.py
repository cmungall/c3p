"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:38390 enone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone of general formula
    R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= H) in which the C=O function
    is conjugated to a C=C double bond at the alpha,beta position.

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

    # Look for enone pattern
    enone_pattern = Chem.MolFromSmarts("[CX3]=[CX3][CX3](=O)[!#1]")
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    
    if not enone_matches:
        return False, "No enone substructure found"

    # Check if the double bond and ketone are conjugated
    for match in enone_matches:
        # Atoms in the match: C=C-C(=O)-R
        c1, c2, c3, r = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Check if the double bond and ketone are in the same pi system
        if not mol.GetAromaticForm().IsBondInRingOfSize(c1.GetBonds()[0], 6):
            continue  # Not conjugated
        
        # Check if R is not hydrogen
        if r.GetAtomicNum() == 1:
            continue  # R is hydrogen
        
        return True, "Contains alpha,beta-unsaturated ketone with conjugated C=C and C=O groups"

    return False, "No valid enone substructure found"