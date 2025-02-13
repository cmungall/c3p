"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define inositol phosphate group with flexible chirality 
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol phosphate group detected"

    # Define ceramide moiety allowing variability in chain attachment
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]([CH2,CH][OX2H1])CO")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide-like moiety detected"
    
    # Phosphodiester linkage between ceramide and inositol
    phosphodiester_pattern = Chem.MolFromSmarts("COP(=O)(O)OCC")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage detected"

    # Long acyl chain characterizing sphingoid base
    long_chain_pattern = Chem.MolFromSmarts("[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long acyl chains detected"

    # Passed all checks, positive identification
    return True, "Contains inositol phosphate group, ceramide moiety, phosphodiester linkage, and appropriate long chains"