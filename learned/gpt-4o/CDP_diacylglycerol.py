"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    A CDP-diacylglycerol is defined as a glycerol molecule with acyl groups
    (unspecified but typically fatty acyl) at the 1- and 2-positions, and a 
    cytidine diphosphate (CDP) group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify glycerol backbone with ester groups at positions 1 and 2
    glycerol_pattern = Chem.MolFromSmarts("[C@H](OC(*)OC(*))COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with esterification at 1,2-positions"

    # Verify acyl groups are unspecified (should typically appear as carbon chains in esters)
    acyl_ester_pattern = Chem.MolFromSmarts("OC(=O)(C*)O[C@@H]")
    acyl_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if len(acyl_matches) < 2:
        return False, "Could not find two unspecified acyl groups (possible esters) at 1,2-positions"
    
    # Identify the CDP group attached at the glycerol's 3-position
    cdp_group_pattern = Chem.MolFromSmarts("n1cc(COP(=O)(O)O[P](=O)(O)*)nc1O")
    if not mol.HasSubstructMatch(cdp_group_pattern):
        return False, "No cytidine diphosphate (CDP) group at glycerol 3-position"
    
    return True, "Contains glycerol backbone with linked acyl at 1,2-positions and CDP group at the 3-position"