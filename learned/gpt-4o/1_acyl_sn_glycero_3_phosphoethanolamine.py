"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphoethanolamine group
    pe_pattern = Chem.MolFromSmarts("O=P(O)(O)OCCN")  # Match exactly with both oxygens on phosphate
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Check for the chiral glycerol backbone
    glycerol_backbone_pattern = Chem.MolFromSmarts("C[C@@H](O)CO")  # Focus on correct stereocenter
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No chiral glycerol backbone found"
    
    # Check for the acyl chain properly linked via ester bond
    acyl_linkage_pattern = Chem.MolFromSmarts("C(=O)OC[C@H]CO")  # Ensure ester link
    if not mol.HasSubstructMatch(acyl_linkage_pattern):
        return False, "No proper acyl linkage to glycerol found"

    # Additional checks can include the long acyl chain characteristics typically found
    long_acyl_chain_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCCCCC")  # Example of a long chain check
    if not mol.HasSubstructMatch(long_acyl_chain_pattern):
        return False, "The acyl chain is not long enough"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoethanolamine structure"