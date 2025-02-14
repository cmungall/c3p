"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    A PE is a glycerol backbone with two fatty acid chains and a phosphate group,
    which is further linked to an ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PE, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone (C-C-C with 3 oxygens attached, two as esters, and one to the phosphate).
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Check for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
         return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for the phosphate group (P(=O)(O)-O)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([OX1])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 1"

    # Check for ethanolamine head group (O-C-C-N)
    ethanolamine_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][NX3]")
    ethanolamine_matches = mol.GetSubstructMatches(ethanolamine_pattern)

    if len(ethanolamine_matches) != 1:
        return False, "Missing ethanolamine head group"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Check number of rotatable bonds - they should be higher than 8
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - PEs typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for PE"

    # Count the number of phosphorus and nitrogen atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"
    if n_count < 1:
       return False, "Must have at least 1 nitrogen atom"
    if o_count < 7:
       return False, "Must have at least 7 oxygen atoms"

    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphate group linked to an ethanolamine"