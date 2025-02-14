"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA(4-) consists of a fatty acid chain with a 3-hydroxy group,
    linked to coenzyme A via a thioester, with 4 negative charges on the phosphate groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core CoA Pattern (without charges)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA core structure found."

    # Check for the thioester link
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester link found."

    # Check for the 3-hydroxy group (simplified pattern)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][CX4](O)[CX4]") # a simpler pattern
    if not mol.HasSubstructMatch(hydroxy_pattern):
         return False, "No 3-hydroxy group found."
    
    # Check for four negatively charged phosphate groups
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 3:
        return False, "Must have 3 phosphorus atoms"
    
    total_negative_charge = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Check if the atom is phosphorus
           for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                    total_negative_charge +=1
    
    if total_negative_charge != 4:
        return False, f"Found {total_negative_charge} negatively charged oxygens on phosphates, need exactly 4."

    # Check for one sulfur atom
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count != 1:
        return False, "Must have 1 sulfur atom"

    # Additional check:
    # Counting the specific atoms - check they are within the expected ranges.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
       return False, "Too few carbon atoms to be a fatty acyl CoA"

    return True, "Molecule matches the criteria for a 3-hydroxy fatty acyl-CoA(4-)."