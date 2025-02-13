"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:35195 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically a long-chain ester formed between a fatty acid and a fatty alcohol,
    characterized by long hydrocarbon chains and minimal other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Must have exactly 1 ester group, found {len(ester_matches)}"

    # Count atoms and check composition
    total_atoms = mol.GetNumAtoms()
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    other_count = total_atoms - c_count - o_count
    
    # Waxes should be mostly carbon and hydrogen
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a wax"
    if o_count != 2:
        return False, f"Waxes should have exactly 2 oxygens (ester group), found {o_count}"
    if other_count > 0:
        return False, "Contains atoms other than C, H, and O"

    # Check for cyclic structures (waxes are typically linear)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings - waxes are typically linear"

    # Look for long carbon chains (at least 8 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Must have at least two long hydrocarbon chains"

    # Count rotatable bonds to verify flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Not enough rotatable bonds for a wax"

    # Check molecular weight - waxes typically >350 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a wax"
    if mol_wt > 1000:
        return False, "Molecular weight too high for typical wax"

    # Check for carboxylic acids
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains free carboxylic acid group"

    # Calculate fraction of sp3 carbons (should be high for waxes)
    sp3_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_count = len(mol.GetSubstructMatches(sp3_pattern))
    sp3_fraction = sp3_count / c_count
    if sp3_fraction < 0.8:
        return False, "Too many unsaturated carbons for a typical wax"

    return True, "Long-chain ester with appropriate molecular weight and composition typical of waxes"