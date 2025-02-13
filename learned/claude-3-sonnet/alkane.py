"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: CHEBI:33677 alkane
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic hydrocarbon with only single bonds (saturated).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains only carbon and hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:  # 1=H, 6=C
            return False, f"Contains non C/H atom: {atom.GetSymbol()}"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings - alkanes must be acyclic"
        
    # Check for unsaturation (double/triple bonds)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:  # 1.0 = single bond
            return False, "Contains unsaturated bonds"
            
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    expected_h = 2*c_count + 2  # CnH2n+2 formula
    
    # Get total hydrogen count (including implicit hydrogens)
    total_h = 0
    for atom in mol.GetAtoms():
        total_h += atom.GetTotalNumHs() + int(atom.GetNumExplicitHs())
    
    # Verify hydrogen count matches formula
    if total_h != expected_h:
        return False, f"H count {total_h} doesn't match expected {expected_h} for alkane formula CnH2n+2"
    
    # Check all carbons are sp3 (saturated)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            if atom.GetHybridization() != Chem.HybridizationType.SP3:
                return False, "Contains non-sp3 carbons"
            
    return True, f"Acyclic saturated hydrocarbon with formula C{c_count}H{total_h}"