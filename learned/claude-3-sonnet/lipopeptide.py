"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI:51842 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for peptide backbone (-C(=O)-N-)
    peptide_pattern = Chem.MolFromSmarts("[C&D2](-[N&D2])=O")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if not peptide_matches:
        return False, "No peptide backbone found"
    
    # Look for lipid chains (long carbon chains)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid chains found"
    
    # Check for attachment of lipid to peptide
    for lipid_match in lipid_matches:
        for peptide_match in peptide_matches:
            if any(atom.GetIdx() in lipid_match for atom in mol.GetAtomWithIdx(peptide_match[0]).GetNeighbors()):
                return True, "Contains a peptide with attached lipid chain"
    
    return False, "No evidence of lipid attached to peptide"