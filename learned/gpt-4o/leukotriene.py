"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is characterized by a C20 polyunsaturated fatty acid backbone with four double bonds, 
    three of which are conjugated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carbon chain typical of icosanoids (around 20 carbons)
    atom_count = mol.GetNumAtoms()
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 20:
        return False, "Molecule does not have adequate icosanoid carbon chain length"
    
    # Count the number of double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches < 4:
        return False, f"Found {double_bond_matches} double bonds, need at least 4"
    
    # Check for conjugated double bonds (at least three consecutive double bonds)
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    conjugated_matches = mol.HasSubstructMatch(conjugated_pattern)
    if not conjugated_matches:
        return False, "No conjugated double bonds pattern found"
    
    # Detect other characteristic groups of leukotrienes (e.g., hydroxyl groups, sulfidopeptide)
    # Let's check for the presence of hydroxyl (OH) as a basic feature
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    num_oh_groups = len(oh_matches)
    
    # Return True with reasons for being classified as leukotriene
    classification_reasons = [
        "Has adequate carbon chain length",
        f"Has {double_bond_matches} double bonds, of which pattern supports conjugation" if conjugated_matches else "",
        f"Contains {num_oh_groups} hydroxyl groups" if num_oh_groups > 0 else "No hydroxyl groups detected"
    ]
    reason = "; ".join(filter(None, classification_reasons))
    
    return True, reason

# Example usage:
# is_leukotriene("CCCCCCCCCCCCCCCCCCCC=O")