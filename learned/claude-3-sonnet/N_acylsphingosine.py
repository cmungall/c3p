"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:38148 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    N-acylsphingosines are composed of sphingosine with an unspecified fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sphingosine backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CHX4][CHX4][CHX3]=[CHX2][CHX2][CHX2][CHX2][CHX2][CHX2][CHX2][CHX2][CHX3]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Look for amide group (-N-C(=O)-) attached to sphingosine
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    amide_match = mol.GetSubstructMatches(amide_pattern)
    if not amide_match:
        return False, "No amide group found"
    
    # Check for fatty acid chain attached to amide
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern, atomDepth=4)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"
    
    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acid chain too short"
    
    # Check molecular weight - N-acylsphingosines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for N-acylsphingosine"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for N-acylsphingosine"
    if o_count != 4:
        return False, "Must have exactly 4 oxygens"
    
    return True, "Contains sphingosine backbone with fatty acyl group attached to nitrogen"