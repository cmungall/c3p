"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with a phosphate group attached via an oxygen
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2][CX3](=[OX1]))[CH2X4][OX2][P](=[OX1])([OX2])")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
         return False, "No glycerol-phosphate backbone found"


    # Check for a myo-inositol ring (six membered ring with 5 OH groups)
    inositol_pattern = Chem.MolFromSmarts("C1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring found"

    # Check for a phosphate group on the inositol ring, directly linked via an oxygen
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C][O][P](=[OX1])([OX2])([OX2])")
    inositol_phosphate_matches = mol.GetSubstructMatches(inositol_phosphate_pattern)
    if len(inositol_phosphate_matches) < 2:  # One phosphate from the glycerol and one additional at least.
        return False, "Less than two phosphate groups found in the phosphoinositide"

    # Check for linkage between glycerol phosphate and inositol phosphate
    phosphodiester_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])[OX2][C]1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester bond found connecting glycerol phosphate to inositol"

    return True, "Contains glycerol backbone, a phosphate group linked to glycerol, a myo-inositol ring, and at least one phosphate group on the inositol ring"