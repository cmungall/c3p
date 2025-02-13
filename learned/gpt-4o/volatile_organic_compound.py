"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is any organic compound with a boiling point <= 250°C under standard atmospheric pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Check the length of the carbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for specific functional groups that significantly affect volatility
    # More relaxed carbon count filter, adjusted by molecular weight for longer chains
    alcohol_pattern = Chem.MolFromSmarts("O")
    ether_pattern = Chem.MolFromSmarts("O-C")
    alkene_pattern = Chem.MolFromSmarts("C=C")
    alkyne_pattern = Chem.MolFromSmarts("C#C")
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]")

    # If there's at least one known VOC functional group and the molecular weight isn't too high
    if (mol.HasSubstructMatch(alcohol_pattern) or 
        mol.HasSubstructMatch(ether_pattern) or 
        mol.HasSubstructMatch(alkene_pattern) or 
        mol.HasSubstructMatch(alkyne_pattern) or 
        mol.HasSubstructMatch(halogen_pattern)) and mol_wt <= 350:
        return True, "Contains VOC functional group and has acceptable molecular weight"

    # Apply different carbon chain logic depending on structure
    if carbon_count <= 15:
        return True, "Carbon chain length and structure suggest volatility"

    return False, "Structure does not indicate VOC features"