"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:33709 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has L-configuration at the alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha-amino acid pattern: N[C@H](*)C(=O)O
    amino_acid_pattern = Chem.MolFromSmarts("N[C@H](*)[C@](=O)O")
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if not matches:
        return False, "No alpha-amino acid substructure found"
    
    # Check for L-configuration at alpha-carbon
    alpha_carbon_idx = matches[0][1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    
    if alpha_carbon.GetProp("_ChiralityPossible") == "No":
        return False, "Alpha-carbon is not a chiral center"
    
    if alpha_carbon.GetProp("_CIPCode") == "S":
        return True, "Contains an L-alpha-amino acid substructure"
    else:
        return False, "Contains a D-alpha-amino acid substructure"