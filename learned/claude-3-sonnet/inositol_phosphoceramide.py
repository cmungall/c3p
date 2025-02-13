"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol ring (cyclohexane with 6 OH groups)
    # More flexible pattern that matches myo-inositol core
    inositol_pattern = Chem.MolFromSmarts("[OX2][CH]1[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]1[OX2]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphodiester bridge
    # Pattern matches P(=O)(O)(OR)(OR) where R can be C or H
    phosphodiester_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])[OX2]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphodiester_pattern))
    if phosphate_matches == 0:
        return False, "No phosphodiester bridge found"
    if phosphate_matches > 1:
        return False, "Multiple phosphate groups found"

    # Check for ceramide core structure (more flexible pattern)
    # Matches both regular and hydroxylated ceramides
    ceramide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide core structure found"

    # Check for long alkyl chains characteristic of ceramides
    alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons in chain
    alkyl_chains = len(mol.GetSubstructMatches(alkyl_chain_pattern))
    if alkyl_chains < 2:
        return False, "Missing long alkyl chains characteristic of ceramides"

    # Check for proper connectivity between components
    # More flexible pattern that captures the essential P-O linkages
    connectivity_pattern = Chem.MolFromSmarts("[CH1]1([OX2][P](=[OX1])([OX2])[OX2]CC[NH])[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]1[OX2]")
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Missing required connectivity between inositol and ceramide"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    # Verify basic composition requirements
    if c_count < 20:
        return False, "Carbon count too low for ceramide chains"
    if o_count < 8:  # 6 from inositol, 1 from amide, minimum 1 from phosphate
        return False, "Insufficient oxygen atoms"
    if n_count != 1:
        return False, "Must contain exactly one nitrogen (ceramide)"
    if p_count != 1:
        return False, "Must contain exactly one phosphorus"

    # Check for optional mannose residue
    mannose_pattern = Chem.MolFromSmarts("O[CH]1O[CH](CO)[CH](O)[CH](O)[CH]1O")
    has_mannose = mol.HasSubstructMatch(mannose_pattern)

    base_description = "Contains inositol ring with phosphodiester-linked ceramide"
    if has_mannose:
        return True, f"{base_description} and mannose residue"
    return True, base_description