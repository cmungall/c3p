"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a tetrahydrochromenochromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core rotenoid skeleton using SMARTS - this will be a simplified core, to catch many substitutions
    # The SMARTS pattern represents the fused ring system with the key oxygen atoms and carbonyl
    # This is a challenging part to get right for all rotenoids.
    core_pattern = Chem.MolFromSmarts("c1cc2OC3c4c(c1)C(=O)C5=C(C=C(O3)C=C5)OC4")
    
    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core tetrahydrochromenochromene skeleton not found"

    # Additional checks could be added here, like for a carbonyl group in a specific position (position 12),
    # number of methoxy groups, if available
    # These additional checks are omitted for simplicity and to classify as many molecules as possible


    return True, "Contains the core tetrahydrochromenochromene skeleton, consistent with a rotenoid"