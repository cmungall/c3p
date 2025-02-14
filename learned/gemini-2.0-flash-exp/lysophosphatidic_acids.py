"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a glycerol backbone with one fatty acid chain attached via an ester bond and a phosphate group.

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

    # Look for glycerol backbone with phosphate pattern (C-C-C with 1 oxygen, 1 phosphate and 1 hydroxyl)
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])([OX2]P)")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for a fatty acid chain (long carbon chain attached to the ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
         return False, f"No fatty acid chain found"

    # Check the number of phosphorus atoms - must be 1.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
         return False, f"Must have exactly 1 phosphorus atom, found {p_count}"
    
    # Count rotatable bonds to verify long chains.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
       return False, "Chain too short to be a fatty acid"
    
    # Check molecular weight - LPA typically > 300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a lysophosphatidic acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 8:
         return False, "Too few carbons for a lysophosphatidic acid"
    if o_count < 6:
         return False, "Must have at least 6 oxygens"

    return True, "Contains a glycerol backbone with a single fatty acid chain and a phosphate group"