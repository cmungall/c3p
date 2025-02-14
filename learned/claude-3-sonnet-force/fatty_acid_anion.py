"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:51710 Fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX1-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"
    
    # Look for long carbon chain
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"
    
    # Check for common functional groups that would disqualify it as a fatty acid anion
    disqualifying_patterns = ["[NX3]", "[#6]#[#7]", "[#6]=[#16]", "[#6]#[#6]"]
    for pattern in disqualifying_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Found disqualifying functional group: {pattern}"
    
    return True, "Contains a carboxylate group and a long carbon chain"