"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: CHEBI:28026 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide is a phosphosphingolipid with an inositol residue and a ceramide
    moiety linked via a phosphodiester bridge.

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

    # Look for inositol moiety
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]([C@@H]([C@@H]([C@H](O)[C@H]1O)O)O)O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol moiety found"

    # Look for phosphodiester bridge (-O-P(=O)(O)-O-)
    phosphodiester_pattern = Chem.MolFromSmarts("[OX2]P(=[OX1])(O)[OX2]")
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if not phosphodiester_matches:
        return False, "No phosphodiester bridge found"

    # Look for ceramide moiety (long fatty acid chain + sphingoid base)
    ceramide_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4]~[NX3][CX3](=[OX1])")
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if not ceramide_matches:
        return False, "No ceramide moiety found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short for ceramide moiety"

    # Check molecular weight - inositol phosphoceramides typically >800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for inositol phosphoceramide"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for inositol phosphoceramide"
    if o_count < 10:
        return False, "Too few oxygens for inositol phosphoceramide"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (sphingoid base)"

    return True, "Contains inositol moiety linked to ceramide via phosphodiester bridge"