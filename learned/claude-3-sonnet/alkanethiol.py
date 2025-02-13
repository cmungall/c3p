"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: A compound in which a sulfanyl group, -SH, is attached to an alkyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for -SH (sulfanyl) group pattern
    sulfanyl_pattern = Chem.MolFromSmarts("[SH1]")
    if not mol.HasSubstructMatch(sulfanyl_pattern):
        return False, "No sulfanyl (-SH) group found"

    # Count sulfanyl groups
    sulfanyl_matches = len(mol.GetSubstructMatches(sulfanyl_pattern))
    if sulfanyl_matches > 2:
        return False, "Too many sulfanyl groups for simple alkanethiol"

    # Exclude aromatic compounds
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Exclude compounds with complex functional groups
    exclusion_patterns = [
        ("[CX3](=O)[OX2H1]", "Contains carboxylic acid"),
        ("[CX3](=O)[OX2][CX4]", "Contains ester group"),
        ("[CX3]=O", "Contains ketone/aldehyde"),
        ("[NX3](~[!#1])(~[!#1])~[!#1]", "Contains tertiary amine"),
        ("[#7X2]=,:[#6]", "Contains imine or similar"),
        ("[SX2](!@[#6])!@[#6]", "Contains thioether bridge"),
        ("[S;!$([SH1]);!$(S(=O));!$(S(=O)=O)]", "Contains other sulfur groups"),
        ("[CX3]=[CX3]", "Contains alkene not adjacent to SH"),
        ("[#6X3]~[#6X3]", "Contains conjugated system"),
        ("[S;R]", "Contains sulfur in ring"),
        ("[C;R5]", "Contains 5-membered ring"),
        ("[C;R6]", "Contains 6-membered ring"),
        ("[PX4]", "Contains phosphorus"),
        ("[BX3]", "Contains boron"),
        ("[SiX4]", "Contains silicon"),
        ("[S;X3,X4,X5,X6]", "Contains oxidized sulfur")
    ]

    for pattern, reason in exclusion_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            return False, reason

    # Get all sulfur atoms and check their environment
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            if atom.GetDegree() == 1:  # Terminal sulfur (as in -SH)
                carbon = atom.GetNeighbors()[0]
                if carbon.GetAtomicNum() != 6:  # Must be connected to carbon
                    return False, "Sulfanyl not connected to carbon"
                if carbon.GetIsAromatic():
                    return False, "Connected to aromatic carbon"
                if carbon.GetHybridization() != Chem.HybridizationType.SP3:
                    return False, "Connected carbon must be sp3 hybridized"

    # Allow simple primary amines (max 1) if present with thiol
    amine_pattern = Chem.MolFromSmarts("[NX3H2]")
    if mol.HasSubstructMatch(amine_pattern):
        if len(mol.GetSubstructMatches(amine_pattern)) > 1:
            return False, "Contains multiple primary amines"
        # Check distance between amine and thiol
        if rdMolDescriptors.CalcNumAtomStereoCenters(mol) > 2:
            return False, "Too complex for simple aminothiol"

    # Final check for molecular complexity
    if rdMolDescriptors.CalcNumRotatableBonds(mol) > 12:
        return False, "Too many rotatable bonds for simple alkanethiol"
    
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings"

    return True, "Contains sulfanyl group (-SH) attached to an alkyl group"