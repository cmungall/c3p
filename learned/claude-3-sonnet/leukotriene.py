"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:36114 leukotriene
A leukotriene is any icosanoid from the family of C20 polyunsaturated fatty acids and their derivatives generated by leukocytes from arachidonic acid, each member having four double bonds of which three are conjugated.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Found {c_count} carbons, leukotriene must have exactly 20"
    
    # Check for 3 conjugated double bonds
    conj_bonds = rdMolDescriptors.CalcNumAromaticRings(mol, 6)
    if conj_bonds != 3:
        return False, f"Found {conj_bonds} conjugated double bonds, leukotriene must have exactly 3"
    
    # Check for carboxylic acid at specific position
    cooh_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCC(=O)O")
    cooh_matches = mol.GetSubstructMatches(cooh_pattern)
    if len(cooh_matches) != 1:
        return False, "Carboxylic acid not found at the expected position"
    
    # Check for alcohol at specific position
    oh_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCC(O)C")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) != 1:
        return False, "Alcohol not found at the expected position"
    
    # Check for common leukotriene substructures like cysteinyl group (optional)
    cys_pattern = Chem.MolFromSmarts("C(N)CS")
    cys_matches = mol.GetSubstructMatches(cys_pattern)
    
    # Check molecular weight range (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight outside the typical range for leukotrienes"
    
    return True, "Molecule meets criteria for leukotriene"