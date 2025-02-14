"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:26119 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is a long-chain organic compound or mixture that is malleable at ambient temperatures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count rotatable bonds as a proxy for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Not enough rotatable bonds, likely too short chain"
    
    # Look for common wax functional groups (esters, ethers, alcohols)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    
    if not (has_ester or has_ether or has_alcohol):
        return False, "No common wax functional group found"
    
    # Check for long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    if not chain_matches:
        return False, "No long carbon chains found"
    
    # Waxes typically have molecular weight > 300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for wax"
    
    return True, "Contains long hydrocarbon chains with ester, ether, or alcohol functional groups"