"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin consists of two phosphatidic acid molecules linked to a glycerol backbone.
    Each phosphatidic acid has two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count phosphorus atoms - must have exactly 2 for cardiolipin
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 2:
        return False, f"Must have exactly 2 phosphorus atoms, found {phosphorus_count}"

    # Count ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 4:
        return False, f"Found {ester_matches} ester groups, need exactly 4 for fatty acid chains"

    # Look for phosphate groups with specific connectivity
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if len(mol.GetSubstructMatches(phosphate_pattern)) != 2:
        return False, "Missing required phosphate groups"

    # Check for glycerol-phosphate connections
    glycerol_phosphate = Chem.MolFromSmarts("[CH2X4]OP(=O)([OX2])[OX2]")
    if len(mol.GetSubstructMatches(glycerol_phosphate)) < 2:
        return False, "Missing glycerol-phosphate connections"

    # Look for fatty acid chains (long carbon chains)
    fatty_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_chain_count = len(mol.GetSubstructMatches(fatty_chain))
    if fatty_chain_count < 4:
        return False, f"Found only {fatty_chain_count} long carbon chains, need at least 4"

    # Count oxygens - cardiolipins should have at least 13 oxygens 
    # (4 esters + 2 phosphates + 3 glycerol backbones)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 13:
        return False, f"Too few oxygens ({oxygen_count}) for cardiolipin structure"

    # Check molecular weight - cardiolipins typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for cardiolipin"

    # Count carbons - cardiolipins typically have >60 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 60:
        return False, f"Too few carbons ({carbon_count}) for typical cardiolipin"

    # Look for central glycerol backbone with two phosphate connections
    central_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]C[CH][CH2][OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(central_pattern):
        # Try alternative pattern with different bond configuration
        alt_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]C[CH2][CH][OX2]P(=O)([OX2])[OX2]")
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "Missing central glycerol-bisphosphate structure"

    return True, "Contains required cardiolipin structure with central glycerol, two phosphates, and four fatty acid chains"