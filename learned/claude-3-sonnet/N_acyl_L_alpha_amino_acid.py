"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for alpha-amino acid core with N-acyl
    # [C:1] is alpha carbon, must have H
    # [N:2] is amino group, must have acyl substituent
    # [C:3] is carbonyl carbon of acyl group
    amino_acid_pattern = Chem.MolFromSmarts("[C:1]([H])([N:2]C(=[O:3])[#6,#7,#8])[C](=[O])[O]")
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if not matches:
        return False, "No N-acyl alpha-amino acid core structure found"
        
    # Get the alpha carbon and check its chirality
    for match in matches:
        alpha_carbon_idx = match[0]
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        
        # Check if carbon is chiral
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            return False, "Unspecified stereochemistry at alpha carbon"
            
        # In RDKit, counterclockwise (CCW) is typically R configuration
        # For amino acids, L configuration corresponds to S configuration at the alpha carbon
        # Due to CIP priority rules, this usually corresponds to CCW/R in RDKit
        if alpha_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "Not L configuration at alpha carbon"
            
        # Check N-acyl group
        nitrogen_idx = match[1]
        carbonyl_idx = match[2]
        
        # Verify N-acyl group
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        if nitrogen.GetTotalNumHs() > 1:  # Should only have 0 or 1 H (not free NH2)
            return False, "No N-acyl group found"
            
    return True, "Contains N-acyl-L-alpha-amino acid structure"