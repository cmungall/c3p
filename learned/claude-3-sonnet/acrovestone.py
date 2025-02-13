"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone compound based on its SMILES string.
    Acrovestone compounds are typically isoflavone derivatives with glycosidic attachments.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chromone core (benzopyran with ketone)
    chromone_pattern = Chem.MolFromSmarts("O=C1C=COc2ccccc12")
    if not mol.HasSubstructMatch(chromone_pattern):
        return False, "No chromone core structure found"

    # Check for presence of sugar moiety patterns
    sugar_pattern = Chem.MolFromSmarts("[OX2][CH]1[CH][CH][CH]([CH][CH]1)(O)[CH2][OX2]")
    if not mol.HasSubstructMatch(sugar_pattern):
        sugar_pattern2 = Chem.MolFromSmarts("[OX2][CH]1[CH][CH][CH]([CH]1)(O)[CH2][OX2]")
        if not mol.HasSubstructMatch(sugar_pattern2):
            return False, "No glycosidic moiety found"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check for aromatic rings
    aromatic_rings = len(rdMolDescriptors.CalcAromaticRings(mol))
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings"

    # Calculate molecular weight - should be substantial due to sugar moieties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for acrovestone compound"

    # Count oxygen atoms - should have multiple due to glycosidic bonds and hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Insufficient oxygen atoms"

    # Look for potential methoxy groups (common in these compounds)
    methoxy_pattern = Chem.MolFromSmarts("COc")
    methoxy_matches = len(mol.GetSubstructMatches(methoxy_pattern))

    # Build reason string
    features = []
    features.append("Contains chromone core")
    features.append(f"Has {hydroxyl_matches} hydroxyl groups")
    features.append(f"Contains {aromatic_rings} aromatic rings")
    features.append("Contains glycosidic moiety")
    if methoxy_matches > 0:
        features.append(f"Has {methoxy_matches} methoxy groups")
    
    return True, "; ".join(features)