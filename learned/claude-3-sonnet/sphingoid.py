"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI:26068 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are defined as 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sphinganine backbone pattern
    sphinganine_pattern = Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH1X4]([CH2X4][CH2X4][CX3](=[OX1])[NX3+,NX3][CH2X4][OX2H,OX1-])[CH1X4]([OX2H,OX1H,OX1-])[CH1X4]([OX2H,OX1H,OX1-])[CH2X4][CH2X4][CH2X4][CH2X4][CH3X4]")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Check for hydroxy and unsaturated derivatives
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H,OX1H,OX1-]")
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(hydroxy_pattern) and not mol.HasSubstructMatch(unsaturated_pattern):
        return False, "Neither hydroxy nor unsaturated derivative"
    
    # Check for long aliphatic chains
    alkyl_chain_pattern = Chem.MolFromSmarts("[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_chain_matches) < 2:
        return False, "Missing long aliphatic chains"
    
    # Check for stereochemistry
    chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
    if not chiral_centers:
        return False, "No stereochemistry specified"
    
    return True, "Contains sphinganine backbone with hydroxy/unsaturated derivatives and long aliphatic chains"