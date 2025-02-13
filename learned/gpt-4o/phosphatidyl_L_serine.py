"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine is a glycerol backbone with two fatty acid chains and a phosphatidyl group
    esterified to the hydroxy group of serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for glycerol backbone with ester links
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester links found"
    
    # Pattern for phosphate group connected serine 
    phosphatidyl_serine_pattern = Chem.MolFromSmarts("P(=O)(O)OCC([N+])C(=O)O")
    if not mol.HasSubstructMatch(phosphatidyl_serine_pattern):
        return False, "No phosphatidyl-L-serine group found"

    # Fatty acid ester link detection
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC")
    if mol.GetSubstructMatches(fatty_acid_pattern) < 2:
        return False, "Less than two fatty acid chains attached"
    
    return True, "Contains glycerol backbone with two fatty acid chains and a phosphatidyl group esterified to serine"