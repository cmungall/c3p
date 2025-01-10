"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax typically contains long alkyl chains connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
   
    # Check for ester group presence with more robust pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkages found"
    
    # Identify long alkyl chains: Recognize both linear and some branched chains
    long_alkyl_pattern = Chem.MolFromSmarts("CCCCCCCCCCC")  # Minimum length of 11 carbons
    if not mol.HasSubstructMatch(long_alkyl_pattern):
        return False, "No sufficient long alkyl chains found"
    
    # Check for molecular weight indicative of larger waxes
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Adjusted threshold to include lighter valid waxes
        return False, "Molecular weight too low for a typical wax compound"
    
    # Consider other structural features typical of waxes
    # For refinement: No complex functional groups are typical of waxes
    non_wax_groups = Chem.MolFromSmarts("[#7,#15,#16,#17,#35,#53]")  # Exclude N, P, S, Cl, Br, I
    if mol.HasSubstructMatch(non_wax_groups):
        return False, "Contains elements atypical of waxes"

    return True, "Contains characteristic ester linkage and long alkyl chains typical of waxes"