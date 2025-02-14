"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    This class is characterized by the formal condensation of the thiol group of coenzyme A
    with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester linkage (C(=O)SCC) - this is part of typical acyl-CoA structures
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found indicative of acyl-CoA"

    # Identify the coenzyme A moiety
    # Try a pattern including parts of the CoA - e.g., nucleoside and phosphate groups
    coa_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing segments indicative of coenzyme A"

    # Look for branching in the fatty acyl chain
    # Check for a carbon with more than 2 bonds that is not at the ends of a chain
    branching_pattern = Chem.MolFromSmarts("[CH2,CH](C)(C)")
    if not mol.HasSubstructMatch(branching_pattern):
        return False, "No branching in the fatty acid chain found"

    return True, "Structure matches branched-chain fatty acyl-CoA"