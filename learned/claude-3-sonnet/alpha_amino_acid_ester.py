"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:35701 alpha-amino acid ester
The amino acid ester derivative obtained the formal condensation of an alpha-amino acid with an alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for alpha-amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)(C(=O)O)[CH2]")
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not matches:
        return False, "No alpha-amino acid backbone found"
    
    # Check for ester group (-O-C(=O)R)
    ester_pattern = Chem.MolFromSmarts("O=C(O[Rb])R")  # [Rb] matches alkyl chains
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check if ester is attached to the alpha-amino acid
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[1])
        for bond in ester_atom.GetBonds():
            if bond.GetOtherAtomIdx(ester_atom.GetIdx()) in matches[0]:
                return True, "Contains alpha-amino acid backbone with ester group attached"
    
    return False, "Ester group not attached to alpha-amino acid"