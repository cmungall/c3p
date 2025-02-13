"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: ChEBI:26829 isoflavone

An isoflavone is defined as any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton 
and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavone(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoflavone skeleton (3-aryl-4H-chromen-4-one)
    isoflavone_pattern = Chem.MolFromSmarts("c1coc2c(c1=O)c(ccc2)-c1ccccc1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Missing isoflavone skeleton"
    
    # Check for aromatic rings
    if not mol.GetAromaticRings():
        return False, "No aromatic rings found"
    
    # Check for carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if len(mol.GetSubstructMatches(carbonyl_pattern)) != 1:
        return False, "Incorrect number of carbonyl groups"
    
    # Check for heterocyclic oxygen
    heterocyclic_oxygen_pattern = Chem.MolFromSmarts("O1ccccc1")
    if len(mol.GetSubstructMatches(heterocyclic_oxygen_pattern)) != 1:
        return False, "Incorrect number of heterocyclic oxygens"
    
    # Check for aryl substituent
    aryl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if len(mol.GetSubstructMatches(aryl_pattern)) != 1:
        return False, "Incorrect number of aryl substituents"
    
    return True, "Contains the 3-aryl-4H-chromen-4-one skeleton characteristic of isoflavones"