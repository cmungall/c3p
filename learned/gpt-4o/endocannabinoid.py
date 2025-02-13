"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids have long hydrocarbon chains (saturated/unsaturated) with linkages to small molecules like ethanolamine or glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for long hydrocarbon chain (saturated/unsaturated)
    long_chain_pattern = Chem.MolFromSmarts("[C;!R]=,:[C;!R]")  # Includes saturated chains
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for common endocannabinoid linkages
    amide_linkage_pattern = Chem.MolFromSmarts("NC(=O)")  # Amide linkage 
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)")  # Ester linkage
    ethanolamine_linkage_pattern = Chem.MolFromSmarts("NCCO")  # Ethanolamine
    glycerol_linkage_pattern = Chem.MolFromSmarts("OCC(O)CO")  # Glycerol linkage
    
    has_amide = mol.HasSubstructMatch(amide_linkage_pattern)
    has_ester = mol.HasSubstructMatch(ester_linkage_pattern)
    has_ethanolamine_linkage = mol.HasSubstructMatch(ethanolamine_linkage_pattern)
    has_glycerol_linkage = mol.HasSubstructMatch(glycerol_linkage_pattern)

    if (has_amide or has_ester) and (has_ethanolamine_linkage or has_glycerol_linkage):
        return True, "Contains structural features typical of endocannabinoids"

    return False, "Does not possess key structural features of endocannabinoids"