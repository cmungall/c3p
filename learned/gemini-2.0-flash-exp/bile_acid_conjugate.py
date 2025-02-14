"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates have a steroidal core and are conjugated to a hydrophilic/charged group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define bile acid core pattern (simplified, may need more complex patterns for edge cases)
    # This pattern captures a 4 ring steroidal core with specified stereochemistry.
    bile_acid_core_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@H]([C@H]([C@@H]([C@@H]3[C@H]1CC[C@@]4([C@@H]3CC[C@@H](C4)C)C)C)C)C")
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
         return False, "No bile acid core structure detected."
     
    # Define common conjugation groups using SMARTS patterns

    # Amino acids (glycine, taurine, alanine, serine, etc.)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]([CX4,OX2])[CX3](=[OX1])")

    # Sulfates
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,H]")

    # Glucuronides
    glucuronide_pattern = Chem.MolFromSmarts("C1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)C(=O)O")
    
    # Sugars
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1")

    # Check for the presence of at least one conjugation group
    if not (mol.HasSubstructMatch(amino_acid_pattern) or \
            mol.HasSubstructMatch(sulfate_pattern) or \
            mol.HasSubstructMatch(glucuronide_pattern) or \
            mol.HasSubstructMatch(sugar_pattern)):
        return False, "No common conjugation group detected."
    
    # Count the number of carbons and oxygens as a sanity check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for a bile acid conjugate"
    
    if o_count < 3:
         return False, "Too few oxygens for a bile acid conjugate"


    return True, "Contains a bile acid core and at least one common conjugation group."