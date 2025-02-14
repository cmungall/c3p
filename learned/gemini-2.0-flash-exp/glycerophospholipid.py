"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerol backbone with a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    glycerol_carbon_atoms = [mol.GetAtomWithIdx(i) for i in glycerol_matches[0]]
    terminal_glycerol_carbons = [glycerol_carbon_atoms[0], glycerol_carbon_atoms[2]]

    # Look for phosphate group connected to glycerol through an ester bond or directly connected to oxygen of the glycerol
    # Allow for variations like -O-P(=O)(O)- , -O-P(=O)(O)O, -O-P(=O)(O)[O-] etc.
    phosphate_pattern1 = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])[OX2]") #  -O-P(=O)(-O)-O
    phosphate_pattern2 = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])") #  -O-P(=O)(-O)-
    phosphate_pattern3 = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])[O-]") #  -O-P(=O)(-O)-[O-]
    phosphate_matches1 = mol.GetSubstructMatches(phosphate_pattern1)
    phosphate_matches2 = mol.GetSubstructMatches(phosphate_pattern2)
    phosphate_matches3 = mol.GetSubstructMatches(phosphate_pattern3)


    phosphate_matches = phosphate_matches1 + phosphate_matches2 + phosphate_matches3

    if not phosphate_matches:
        return False, "No phosphate group found"
        
    #Check if the phosphate group is connected to a terminal carbon in glycerol
    phosphate_attached = False
    for match in phosphate_matches:
            phosphorous_atom = mol.GetAtomWithIdx(match[1])
            for terminal_carbon in terminal_glycerol_carbons:
                for neighbor in terminal_carbon.GetNeighbors():
                    if neighbor.GetIdx() == phosphorous_atom.GetNeighbors()[0].GetIdx() and neighbor.GetSymbol() == 'O':
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
            if phosphate_attached:
                break    
    if not phosphate_attached:
        return False, "Phosphate not attached to a terminal glycerol carbon"

    # Look for at least one ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"


   # Check for fatty acid chains attached to glycerol through ester bond
    fatty_acid_ester_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]~[OX2][CX3](=[OX1])") # chain connected to ester
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_ester_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chains found"


    return True, "Contains glycerol backbone with at least one fatty acid and a phosphate group linked via an ester bond to terminal glycerol carbon."