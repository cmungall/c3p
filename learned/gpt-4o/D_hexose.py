"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose (6 C sugar) with D-configuration at the 5th carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has 6 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"
    
    # Check for hydroxyl groups (–OH), more than 3 expected in a hexose
    oh_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))
    if len(oh_groups) < 4:
        return False, f"Expected at least 4 OH groups, found {len(oh_groups)}"
    
    # 5-membered (furanose) or 6-membered (pyranose) ring detection with D-config
    furanose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)C1O")
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C1O")

    # Search for furanose or pyranose pattern and correct stereo config
    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        for atom in mol.GetAtoms():
            # Look for specific D-configuration common at position 5
            if atom.GetChiralTag() and atom.GetAtomicNum() == 6:
                # Find position of D-configuration 
                # Applies to specific stereocenters in common D-hexoses
                chiral_smarts = "[C@H](O)"
                if any(mol.GetAtomWithIdx(idx).HasSubstructMatch(Chem.MolFromSmarts(chiral_smarts)) for idx in range(len(mol.GetAtoms()))):
                    return True, "Identified as D-hexose with correct stereochemistry"
    
    return False, "No matching D-hexose ring pattern found"