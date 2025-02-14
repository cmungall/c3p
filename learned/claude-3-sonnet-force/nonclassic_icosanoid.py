"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:36738 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signalling molecule made by oxygenation of C20 fatty acids
    other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Must have exactly 20 carbon atoms
    if c_count != 20:
        return False, f"Molecule has {c_count} carbon atoms, nonclassic icosanoids must have 20"
    
    # Look for carboxyl group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"
    
    # Look for multiple double bonds
    double_bond_pattern = Chem.MolFromSmarts("=")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, "Fewer than 3 double bonds found"
    
    # Look for at least 3 oxygens
    if o_count < 3:
        return False, "Fewer than 3 oxygen atoms found"
    
    # Check molecular weight - nonclassic icosanoids typically 300-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, "Molecular weight outside typical range for nonclassic icosanoids (300-500 Da)"

    return True, "Contains 20 carbons, carboxyl group, multiple double bonds, and multiple oxygens"