"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Identify carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]

    # Check that all carbons are sp3 hybridized (single bonds)
    for atom in carbons:
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False, "Carbon atom not sp3 hybridized"

    # Identify terminal carbons (carbons connected to only one other carbon)
    terminal_carbons = []
    for atom in carbons:
        num_carbons = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        if num_carbons == 1:
            terminal_carbons.append(atom)

    if len(terminal_carbons) != 2:
        return False, "Molecule does not have exactly two terminal carbons"

    # Check that terminal carbons are -CH2OH groups
    for atom in terminal_carbons:
        neighbors = atom.GetNeighbors()
        has_oh = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if oxygen is part of hydroxyl group (O connected to H)
                if len(neighbor.GetNeighbors()) == 1:
                    has_oh = True
        if not has_oh:
            return False, "Terminal carbon does not have hydroxyl group"

    # Check internal carbons
    for atom in carbons:
        if atom not in terminal_carbons:
            # Should be connected to two carbons and one hydroxyl group
            num_carbons = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if num_carbons != 2:
                return False, "Internal carbon not connected to exactly two carbons"
            has_oh = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    # Check if oxygen is part of hydroxyl group (O connected to H)
                    if len(neighbor.GetNeighbors()) == 1:
                        has_oh = True
            if not has_oh:
                return False, "Internal carbon does not have hydroxyl group"

    # Check for linear chain (no branching)
    for atom in carbons:
        num_carbons = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        if atom in terminal_carbons:
            if num_carbons != 1:
                return False, "Terminal carbon connected to more than one carbon"
        else:
            if num_carbons != 2:
                return False, "Internal carbon connected to incorrect number of carbons"

    return True, "Molecule is an alditol (acyclic polyol with terminal CH2OH groups)"