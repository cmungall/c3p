"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone, a phosphate group, and two fatty acid chains.
    The two fatty acid chains are attached via ester bonds to the glycerol.
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with phosphate attached
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]([OX2][P](=[OX1])([OX2])([OX2]))[CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with attached phosphate found"

    # Look for 2 ester groups attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX3](=[OX1])") # ester to glycerol
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2 :
         return False, f"Found {len(ester_matches)} ester groups attached to glycerol, need exactly 2"
    
    #Check for long fatty acid chains attached to the glycerol via ester
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX3](=[OX1])-[CX4,CX3]~[CX4,CX3]") # fatty acid (>= 2 carbons) attached via ester
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 2:
        return False, f"Missing two fatty acid chains attached via ester bonds, found {len(fatty_acid_matches)}"

    # Look for the minimum required structure
    minimum_pattern = Chem.MolFromSmarts("[CH2X4]([OX2][P](=[OX1])([OX2])([OX2]))[CHX4](-[OX2][CX3](=[OX1])-[CX4])-[CH2X4]-[OX2][CX3](=[OX1])-[CX4]")
    if not mol.HasSubstructMatch(minimum_pattern):
         return False, "Missing minimum structure (glycerol, phosphate, 2 esters)."


    return True, "Contains glycerol backbone with a phosphate group and two fatty acid chains attached via ester bonds"