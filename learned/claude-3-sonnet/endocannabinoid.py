"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:77824 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are cannabinoids present in mammalian biological fluids and tissues that activate cannabinoid receptors.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ethanolamine group (NCCO)
    ethanolamine_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+]CCCO")
    ethanolamine_match = mol.GetSubstructMatches(ethanolamine_pattern)
    
    # Look for glycerol group (C(CO)CO)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    glycerol_match = mol.GetSubstructMatches(glycerol_pattern)
    
    # Look for long aliphatic chain (>= C14)
    aliphatic_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    aliphatic_match = mol.GetSubstructMatches(aliphatic_pattern)
    
    # Look for unsaturated aliphatic chain
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    unsaturated_match = mol.GetSubstructMatches(unsaturated_pattern)
    
    # Check for amide group (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)[NX3]")
    amide_match = mol.GetSubstructMatches(amide_pattern)
    
    # Classify as endocannabinoid if it has:
    # - Ethanolamine or glycerol group
    # - Long aliphatic chain (>= C14)
    # - Unsaturated aliphatic chain
    # - Amide group
    if (ethanolamine_match or glycerol_match) and aliphatic_match and unsaturated_match and amide_match:
        return True, "Contains ethanolamine/glycerol, long aliphatic chain, unsaturation, and amide group"
    else:
        return False, "Missing one or more key structural features of endocannabinoids"