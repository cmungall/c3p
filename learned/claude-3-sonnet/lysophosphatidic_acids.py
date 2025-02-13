"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:75856 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has:
    - A glycerol backbone
    - One fatty acid chain attached via ester bond
    - A phosphate group at the sn-3 position
    - A free hydroxyl group
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for complete lysophosphatidic acid core structure
    # [OH]-C-C(OH/OR)-C-O-P(=O)(O)(O)
    lpa_pattern = Chem.MolFromSmarts(
        "[OX2H,OX2][CH2X4][CHX4]([OX2H,OX2C])[CH2X4][OX2]P(=[OX1])([OX2H,OX1-])[OX2H,OX1-]"
    )
    if not mol.HasSubstructMatch(lpa_pattern):
        return False, "Missing core lysophosphatidic acid structure"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 1:
        return False, f"Found {ester_matches} ester groups, need exactly 1"

    # Look for fatty acid chain attached to ester
    # At least 4 carbons in chain
    fatty_chain = Chem.MolFromSmarts("[CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No suitable fatty acid chain found"

    # Count key atoms to ensure proper composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 5:
        return False, "Too few carbons for lysophosphatidic acid"
    if o_count < 6:
        return False, "Must have at least 6 oxygens"
    if p_count != 1:
        return False, "Must have exactly one phosphorus"

    # Verify at least one free hydroxyl
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "No free hydroxyl group found"

    # Check that molecule isn't too complex (to avoid larger phospholipids)
    if len(mol.GetAtoms()) > 100:
        return False, "Molecule too large for lysophosphatidic acid"

    return True, "Contains glycerol backbone with one fatty acid chain and phosphate group"