"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has an L-carnitine backbone with an acyl group attached via an ester bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic carnitine backbone structure:
    # [C@H] or [C@@H] central carbon with:
    # - CH2 with carboxylate
    # - CH2 with trimethylammonium
    # - Oxygen (part of ester)
    carnitine_pattern = Chem.MolFromSmarts('[C@H,C@@H]([CH2][C]([O-])=O)([CH2][N+](C)(C)C)[OX2]')
    
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "Missing carnitine backbone structure"
        
    # Check for ester linkage
    # The oxygen from carnitine should be connected to a C(=O)R group
    ester_pattern = Chem.MolFromSmarts('[C@H,C@@H]([CH2][C]([O-])=O)([CH2][N+])O[C](=O)')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing or incorrect ester linkage"

    # Check for trimethylammonium group
    # Allow for deuterated methyl groups ([CH3,CD3])
    trimethyl_pattern = Chem.MolFromSmarts('[N+]([CH3,CD3])([CH3,CD3])([CH3,CD3])')
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Missing trimethylammonium group"
    
    # Verify charges
    pos_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    if pos_charge != 1 or neg_charge != -1:
        return False, f"Incorrect charge distribution: +{pos_charge}, {neg_charge}"

    # Check stereochemistry
    # Get all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "Missing required stereocenter"
    
    # Find the central carbon (the one with OH, CH2COO-, and CH2N+)
    central_matches = mol.GetSubstructMatches(carnitine_pattern)
    if not central_matches:
        return False, "Cannot determine stereochemistry"
        
    # The central carbon should have R configuration for L-carnitine
    # (Note: The R configuration gives L-carnitine due to CIP priority rules)
    central_carbon_idx = central_matches[0][0]
    found_correct_config = False
    
    for idx, config in chiral_centers:
        if idx == central_carbon_idx:
            if config == 'R':  # R configuration corresponds to L-carnitine
                found_correct_config = True
            break
            
    if not found_correct_config:
        return False, "Incorrect stereochemistry - must be L-configuration"

    return True, "Valid O-acyl-L-carnitine structure with correct L-configuration"