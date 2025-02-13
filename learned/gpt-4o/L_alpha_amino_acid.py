"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has an amino group, carboxyl group, and an alkyl side chain with the L-configuration at the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to recognize the L-alpha-amino acid backbone
    pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)O")
    
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Structure does not match the L-alpha-amino acid backbone"

    # Check the stereochemistry of the alpha carbon
    for match in matches:
        alpha_carbon = mol.GetAtomWithIdx(match[0])
        chiral_tag = alpha_carbon.GetChiralTag()
        
        # CHI_TETRAHEDRAL_CCW indicates 'S' stereochemistry in RDKit's default settings
        if chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return True, "Contains L-alpha-amino acid structure with L-configuration"
        
    return False, "Structure does not have the L-configuration at the alpha carbon"