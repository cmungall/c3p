"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all rings in the molecule
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "No ring information found in the molecule"

    # Identify sugar rings (5 or 6-membered rings with oxygen)
    pyranose = Chem.MolFromSmarts("OC1CC(O)C(O)C(O)C1O")
    furanose = Chem.MolFromSmarts("OC1CC(O)C(O)C1O")
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1][OR1]1")  # Simplified sugar ring

    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar rings found"

    # Count amino sugars
    amino_sugar_count = 0
    for match in sugar_matches:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Check for amino group attachment (-NH2)
        has_amino_group = False
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                        if neighbor.GetTotalDegree() == 1:
                            has_amino_group = True
                            break
                if has_amino_group:
                    break
        if has_amino_group:
            amino_sugar_count += 1

    if amino_sugar_count == 0:
        return False, "No amino sugars found"
    
    # Check if a substantial proportion of sugars are amino sugars
    total_sugars = len(sugar_matches)
    proportion = amino_sugar_count / total_sugars

    if proportion < 0.5:
        return False, f"Only {amino_sugar_count} out of {total_sugars} sugars are amino sugars"
    
    return True, f"Contains {amino_sugar_count} amino sugars out of {total_sugars} sugar units"