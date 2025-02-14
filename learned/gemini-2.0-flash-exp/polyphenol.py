"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol has 2 or more benzene rings, each with at least one directly attached -OH group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polyphenol, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define benzene and hydroxyl group substructure
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    hydroxyl_pattern = Chem.MolFromSmarts("[c]O")

    # Find all benzene rings
    benzene_rings = mol.GetSubstructMatches(benzene_pattern)
    if not benzene_rings:
         return False, "Molecule does not contain any benzene rings."

    valid_ring_count = 0
    for ring in benzene_rings:
        # Create a query to search for a hydroxyl directly attached to specific atoms of the ring
        benzene_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        ring_has_oh = False
        for atom in benzene_atoms:
          if mol.HasSubstructMatch(hydroxyl_pattern, useChirality=False):
            for match in mol.GetSubstructMatches(hydroxyl_pattern, useChirality=False):
                if mol.GetAtomWithIdx(match[0]).GetNeighbors()[0].GetIdx() == atom.GetIdx():
                  ring_has_oh = True
                  break
          if ring_has_oh:
            break
        if ring_has_oh:
          valid_ring_count += 1


    # Classify based on ring count
    if valid_ring_count >= 2:
        return True, "Contains two or more benzene rings, each with at least one directly attached -OH group."
    else:
        return False, f"Found {len(benzene_rings)} benzene rings, but only {valid_ring_count} have at least one -OH attached, needs at least two."