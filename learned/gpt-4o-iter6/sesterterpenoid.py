"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are derived from sesterterpenes and can have a modified C25 skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 23 or c_count > 50:
        return False, f"Carbon count ({c_count}) not within typical range for sesterterpenoid"

    # Look for complex cyclic structures typical in terpenoids/sesterterpenoids
    complex_cyclic_patterns = [
        Chem.MolFromSmarts("C1CCC2CCC(C1)C2"),  # Simple bicyclic system
        Chem.MolFromSmarts("C1CCCC2CCCC3(C1)C2CCC3"),  # Tricyclic system
    ]

    # Check for presence of complex cyclic structures typical in sesterterpenoids
    if not any(mol.HasSubstructMatch(pattern) for pattern in complex_cyclic_patterns):
        return False, "No key cyclic structures typical of sesterterpenoids found"
    
    # Calculate and validate additional characteristics like molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:  # Adjusting expected range based on known compounds
        return False, "Molecular weight not within typical range for sesterterpenoids"

    # If structure fits broad criteria
    return True, "Likely a sesterterpenoid based on carbon count and presence of key cyclic structures"