"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerol backbone with a phosphate group
    attached to the 3-position, esterified with ethanolamine and a fatty acid at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for the glycerol backbone with phosphate and ethanolamine substitution
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]([OX2][CX3](=[OX1])[#6])[CHX4](O)[CH2X4]OP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with phosphate, ethanolamine and ester at position 1 not found"
    
    # Check for fatty acid chain (long carbon chain) at the 1-position
    # Looking for at least 4 carbons attached to the ester at position 1 - now included in the glycerol pattern
    # Count rotatable bonds (for fatty acid chain check)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acid chain too short"
    
    #Check for at least one phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
       return False, "Must have exactly one phosphorus"

    # Count oxygens - now we expect at least 6, potentially more if the fatty acid contains more
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
       return False, "Must have at least 6 oxygens, includes oxygens in glycerol, ester, phosphate and ethanolamine"


    return True, "Contains glycerol backbone with phosphate and ethanolamine and a fatty acid chain at the 1-position"