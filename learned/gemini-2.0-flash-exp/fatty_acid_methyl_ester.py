"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid with the carboxyl group esterified with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the methyl ester group pattern (C(=O)OC)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"

    # Look for fatty acid chains (long carbon chains attached to carbonyl of ester)
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "Missing fatty acid chain"

    # Verify chain length using rotatable bonds, ensure at least 3 (for a chain of 4)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Fatty acid chain too short"
    
    # Check for glycerol backbone pattern to ensure it is not a triglyceride
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone found. Not a fatty acid methyl ester"


    return True, "Contains a methyl ester group and a fatty acid chain"