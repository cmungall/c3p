"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:27439 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton (pentacyclic triterpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hopane skeleton SMARTS pattern
    hopane_smarts = "C1CCC2(C1)CCC3(C2)CCC4(C3)CCCC5(C4)CCCC5C"
    hopane_pattern = Chem.MolFromSmarts(hopane_smarts)
    if hopane_pattern is None:
        return None, "Invalid SMARTS pattern for hopane skeleton"

    # Check for the hopane skeleton in the molecule
    if not mol.HasSubstructMatch(hopane_pattern):
        return False, "No hopane skeleton found"

    # Additional checks (optional): Verify triterpenoid characteristics
    # Count the number of rings (should be 5 for hopanoids)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings ({num_rings}) for a hopanoid"

    # Check molecular weight range for triterpenoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt} Da) outside typical range for hopanoids"

    return True, "Contains hopane skeleton characteristic of hopanoids"