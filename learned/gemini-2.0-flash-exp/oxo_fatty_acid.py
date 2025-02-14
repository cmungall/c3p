"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid has a carboxylic acid group and at least one additional ketone or aldehyde group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
       return False, "No carboxylic acid group found"
    
    #Check for extra carboxylic acid
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxyl_matches) > 1:
        return False, "More than one carboxylic acid group found"

    # Check for ketone or aldehyde group (excluding the carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)") # Explicit Hydrogen needed

    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    if not (has_ketone or has_aldehyde):
       return False, "No ketone or aldehyde group found"

    # Check for a long carbon chain (at least 4 carbons)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
         return False, "Carbon chain too short"
         
    # Molecular weight check for fatty acid 
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
      return False, "Molecular weight too low"
      
    # Carbon count to check for chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Too few carbons"

    return True, "Contains a carboxylic acid, a ketone or aldehyde, and a long carbon chain"