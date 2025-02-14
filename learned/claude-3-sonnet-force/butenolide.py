"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:30335 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone consisting of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-furanone ring pattern
    furanone_pattern = Chem.MolFromSmarts("[O,X1]=[C,X3]-[C,X3]-[C,X3]-[C,X2]")
    if not mol.HasSubstructMatch(furanone_pattern):
        return False, "No 2-furanone ring found"
    
    # Look for lactone (cyclic ester) functional group
    lactone_pattern = Chem.MolFromSmarts("[O,X2]=[C,X3]-[O,X1]-[C,X3]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Count number of rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings != 1:
        return False, f"Expected 1 ring, found {n_rings}"
    
    # Check that the ring is 5-membered (gammma-lactone)
    ring_info = mol.GetRingInfo()
    if ring_info.AtomRings()[0].Count() != 5:
        return False, "Ring is not a 5-membered gammma-lactone"
    
    # Check for double bond in the ring
    double_bonds = mol.GetBonds()
    has_ring_double_bond = any(bond.GetBondType() == Chem.BondType.DOUBLE and 
                                bond.IsInRing() for bond in double_bonds)
    if not has_ring_double_bond:
        return False, "No double bond in the ring"
    
    # Additional check for molecular formula containing carbon, hydrogen and oxygen only
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    allowed_atoms = set("CHO")
    atoms = set(formula.split())
    if not atoms.issubset(allowed_atoms):
        return False, "Molecule contains atoms other than C, H and O"
    
    return True, "Molecule matches the 2-furanone ring with a lactone ring, fulfilling the criteria for a butenolide"