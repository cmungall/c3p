"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:38685 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are a diverse group of secondary metabolites characteristic to Cannabis plant,
    containing oxygen as part of a heterocyclic ring or in various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for characteristic cannabinoid ring systems (benzopyran, pyran, etc.)
    ring_patterns = [Chem.MolFromSmarts(pat) for pat in ["c1ccc2c(c1)OCO2", "c1ccc2c(c1)OC=CO2", "c1ccc2c(c1)OCC=C2"]]
    ring_matches = [mol.HasSubstructMatch(pat) for pat in ring_patterns]
    if not any(ring_matches):
        return False, "No characteristic cannabinoid ring system found"

    # Look for long aliphatic chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long aliphatic chain found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short for cannabinoid"

    # Check for oxygen-containing functional groups (ester, ether, hydroxyl, etc.)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Not enough oxygen-containing functional groups"

    # Check for aromatic rings
    aromatic_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if ring.IsAromatic()]
    if not aromatic_rings:
        return False, "No aromatic rings found"

    return True, "Contains characteristic cannabinoid ring system, long aliphatic chain, and oxygen-containing functional groups"