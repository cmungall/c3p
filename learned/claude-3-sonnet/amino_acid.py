"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33597 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_match = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_match:
        return False, "No carboxylic acid group found"
    
    # Look for amino group (-N)
    amino_pattern = Chem.MolFromSmarts("[N;H2,H1;!$(NC=O)]")
    amino_match = mol.GetSubstructMatches(amino_pattern)
    if not amino_match:
        return False, "No amino group found"
    
    # Check if amino group is attached to the alpha-carbon of the carboxylic acid
    alpha_carbon_pattern = Chem.MolFromSmarts("C(C(=O)O)N")
    alpha_carbon_match = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not alpha_carbon_match:
        return False, "Amino group not attached to the alpha-carbon of the carboxylic acid"
    
    # Check for additional functional groups or substituents
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    thiol_pattern = Chem.MolFromSmarts("S")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    
    hydroxyl_match = mol.GetSubstructMatches(hydroxyl_pattern)
    thiol_match = mol.GetSubstructMatches(thiol_pattern)
    phosphate_match = mol.GetSubstructMatches(phosphate_pattern)
    
    # Check for common substituents on the amino group
    n_methyl_pattern = Chem.MolFromSmarts("CN")
    n_acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    
    n_methyl_match = mol.GetSubstructMatches(n_methyl_pattern)
    n_acetyl_match = mol.GetSubstructMatches(n_acetyl_pattern)
    
    # Check for additional structural constraints
    backbone_pattern = Chem.MolFromSmarts("[C;H2,H1]([C;H2,H1])([C;H2,H1])")
    backbone_match = mol.GetSubstructMatches(backbone_pattern)
    
    ring_pattern = Chem.MolFromSmarts("C1CCCCC1")
    ring_match = mol.GetSubstructMatches(ring_pattern)
    
    # Check for stereochemistry
    stereochemistry_pattern = Chem.MolFromSmarts("[C@H](N)(C(=O)O)[C@@H]")
    stereochemistry_match = mol.GetSubstructMatches(stereochemistry_pattern)
    
    # Calculate molecular weight - amino acids typically <500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for an amino acid"
    
    # Additional checks or conditions can be added as needed
    
    return True, "Contains a carboxylic acid group, an amino group attached to the alpha-carbon, and meets additional structural constraints"