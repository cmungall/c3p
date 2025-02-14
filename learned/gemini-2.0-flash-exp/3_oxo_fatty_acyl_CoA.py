"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA has a CoA moiety, a 3-oxo group and a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety SMARTS (simplified, focusing on key linkage)
    # The core part is the phospho-pantetheine moiety linked to the adenosine 5'-phosphate part
    # Note that this might still need refinement for better specificity
    coa_smarts = "[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define the 3-oxo group SMARTS
    oxo_smarts = "[CX3](=[OX1])[CH2][CX3](=[OX1])"
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)

    if len(oxo_matches) == 0:
       return False, "No 3-oxo group found"
    elif len(oxo_matches) > 1:
      return False, "More than one 3-oxo group found"

    # Check for thioester group linked to CoA
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    if len(thioester_matches) == 0:
         return False, "No thioester group found"

    # Verify that the thioester is close to the coA part
    # We get the indices of the S and C of the thioester, then we search for the indices of the P in the CoA part
    # If these are far away from each other in the graph, then they might not be close
    for thioester_match in thioester_matches:
       thioester_S_idx = thioester_match[1]
       thioester_C_idx = thioester_match[0]
       
    for coa_match in mol.GetSubstructMatches(coa_pattern):
       coa_P_indices = [idx for idx in coa_match if mol.GetAtomWithIdx(idx).GetSymbol() == 'P']
       found_close = False
       for coa_P_idx in coa_P_indices:
           # Check if any of them is within a given distance
           path = Chem.GetShortestPath(mol, thioester_S_idx, coa_P_idx)
           if path is not None and len(path) < 15:
             found_close = True
             break
       if not found_close:
         return False, "Thioester and CoA not connected"


    # Check for a long carbon chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acyl chain found"
    
    # Count rotatable bonds to verify long chains, this time is more relaxed
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be fatty acids"


    # Check molecular weight - acyl-CoA typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, "Molecular weight too low for acyl-CoA"

    return True, "Contains CoA moiety, 3-oxo group, and a fatty acyl chain"