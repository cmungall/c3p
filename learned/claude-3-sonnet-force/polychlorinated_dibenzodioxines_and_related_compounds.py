"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin or related compound.
    These compounds include:
    - Polychlorinated dibenzodioxins
    - Polychlorinated dibenzofurans
    - Polychlorinated biphenyls
    - Polybrominated biphenyls
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxin or related compound
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define simpler core structure patterns
    dibenzodioxin_pattern = Chem.MolFromSmarts('c1ccc2Oc3ccccc3Oc2c1')  # Basic dibenzodioxin
    dibenzofuran_pattern = Chem.MolFromSmarts('c1ccc2oc3ccccc3c2c1')    # Basic dibenzofuran
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')          # Basic biphenyl
    
    # Count halogens
    halogen_pattern = Chem.MolFromSmarts('[Cl,Br]')
    if halogen_pattern is None:
        return False, "Error in halogen SMARTS pattern"
    
    halogen_matches = mol.GetSubstructMatches(halogen_pattern)
    num_halogens = len(halogen_matches)
    
    if num_halogens == 0:
        return False, "No halogens present"

    # Count chlorines and bromines separately
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))

    # Check for core structures
    core_found = False
    core_type = ""
    
    if dibenzodioxin_pattern is not None and mol.HasSubstructMatch(dibenzodioxin_pattern):
        core_found = True
        core_type = "dibenzodioxin"
    elif dibenzofuran_pattern is not None and mol.HasSubstructMatch(dibenzofuran_pattern):
        core_found = True
        core_type = "dibenzofuran"
    elif biphenyl_pattern is not None and mol.HasSubstructMatch(biphenyl_pattern):
        core_found = True
        core_type = "biphenyl"

    if not core_found:
        return False, "No characteristic core structure found"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings"

    # Check for natural product characteristics
    # Count oxygen-containing substituents (excluding core oxygens)
    oh_pattern = Chem.MolFromSmarts('[OH]')
    ome_pattern = Chem.MolFromSmarts('[OC]')  # Matches any O-C bond
    
    oh_count = len(mol.GetSubstructMatches(oh_pattern)) if oh_pattern else 0
    ome_count = len(mol.GetSubstructMatches(ome_pattern)) if ome_pattern else 0
    
    # Adjust for core oxygens in dioxins/furans
    if core_type in ["dibenzodioxin", "dibenzofuran"]:
        ome_count = max(0, ome_count - (2 if core_type == "dibenzodioxin" else 1))

    # Natural product filter - allow some substitution but not excessive
    if oh_count + ome_count > 6:
        return False, "Too many oxygen substituents - likely a complex natural product"

    # Minimum halogen requirements
    if core_type == "dibenzodioxin" and num_halogens < 4:
        return False, "Insufficient halogenation for dibenzodioxin"
    elif core_type == "dibenzofuran" and num_halogens < 4:
        return False, "Insufficient halogenation for dibenzofuran"
    elif core_type == "biphenyl" and num_halogens < 2:
        return False, "Insufficient halogenation for biphenyl"

    # Build classification reason
    halogen_desc = []
    if num_cl > 0:
        halogen_desc.append(f"{num_cl} chlorine{'s' if num_cl > 1 else ''}")
    if num_br > 0:
        halogen_desc.append(f"{num_br} bromine{'s' if num_br > 1 else ''}")
    halogen_str = " and ".join(halogen_desc)
    
    reason = f"Contains {core_type} core with {halogen_str}"
    
    return True, reason