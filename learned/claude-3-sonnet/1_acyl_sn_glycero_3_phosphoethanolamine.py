"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: CHEBI:18013 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoethanolamine is a glycerol backbone with an acyl chain attached
    at the sn-1 position and a phosphoethanolamine group at the sn-3 position, with (R)-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for acyl chain at sn-1 position (-O-C(=O)-) with long carbon chain
    acyl_chain_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_chain_matches) != 1:
        return False, f"Found {len(acyl_chain_matches)} acyl chains, need exactly 1"

    # Look for phosphoethanolamine group at sn-3 position
    pe_group_pattern = Chem.MolFromSmarts("[OX2]P([OX2][CX4][CX4][NX3])(=[OX1])[OX2]")
    pe_group_matches = mol.GetSubstructMatches(pe_group_pattern)
    if len(pe_group_matches) != 1:
        return False, f"Found {len(pe_group_matches)} phosphoethanolamine groups, need exactly 1"

    # Check for (R)-configuration
    chiral_centers = Chem.FindMolChiralUnspecifiedUnbracketedCenters(mol)
    if len(chiral_centers) != 1:
        return False, "Expected exactly 1 chiral center (glycerol backbone)"
    
    chiral_atom = mol.GetAtomWithIdx(chiral_centers[0])
    if chiral_atom.GetProp("_CIPCode") != "R":
        return False, "Glycerol backbone not in (R)-configuration"

    return True, "Contains glycerol backbone with acyl chain at sn-1 and phosphoethanolamine at sn-3, in (R)-configuration"