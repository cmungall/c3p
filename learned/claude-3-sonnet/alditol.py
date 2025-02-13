"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:28036 alditol
An alditol is a carbohydrate that is an acyclic polyol having the general formula 
HOCH2[CH(OH)]nCH2OH (formally derivable from an aldose by reduction of the carbonyl group).
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for acyclic structure
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, alditols must be acyclic"
    
    # Check for hydroxyl (-OH) and ether (-O-) groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(hydroxyl_matches) < 2 or len(ether_matches) < 1:
        return False, "Insufficient hydroxyl and/or ether groups"
    
    # Check for alternating -CH2-, -CH(OH)- pattern
    alditol_pattern = Chem.MolFromSmarts("[CH2][CH]([OH])[CH2]")
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "Does not match alditol pattern HOCH2[CH(OH)]nCH2OH"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3 or o_count < 3:
        return False, "Too few carbon and/or oxygen atoms for an alditol"
    
    return True, "Matches alditol pattern HOCH2[CH(OH)]nCH2OH"