"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is a glycerol backbone with two fatty acid chains and a
    phosphoglycerol group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphatidylglycerol core structure in one SMARTS pattern
    # This pattern matches the glycerol with two ester groups, a phosphate, and another glycerol.
    # [CH2][CH]([CH2]O[P](=O)(O)-O[CH2][CH](O)[CH2]O)
    pg_core_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4][OX2][P](=[OX1])([OX2])[OX2][CH2X4][CHX4]([OX2])[CH2X4][OX2])")
    if not mol.HasSubstructMatch(pg_core_pattern):
        return False, "Phosphatidylglycerol core not found"

    # Check for two fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(ester_matches) < 2 or len(fatty_acid_matches) < 2:
      return False, f"Must have at least 2 ester groups and 2 fatty acid chains"

    #Check rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 : #Reduced threshold to 8 to accommodate shorter chains
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - PGs are typically >500 Da, but reducing to allow smaller PGs
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for phosphatidylglycerol"
    

    # Count phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 1:
        return False, "Must have at least one phosphorus"

    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphoglycerol group"