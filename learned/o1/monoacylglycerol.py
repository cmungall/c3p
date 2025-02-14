"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which any one of the R groups (position not specified)
    is an acyl group while the remaining two R groups can be either H or alkyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern (three carbons each possibly connected to oxygen)
    glycerol_pattern = Chem.MolFromSmarts("[CX4H2,O][CX4H][CX4H2,O]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Define ester linkage pattern (acyl group attached via ester bond)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Identify acyl group
    acyl_carbon = None
    for match in ester_matches:
        # The carbonyl carbon in the ester group
        acyl_carbon = match[0]
        break

    # Check that the ester group is connected to the glycerol backbone
    glycerol_atoms = set(sum(glycerol_matches, ()))
    if acyl_carbon is not None:
        paths = Chem.rdmolops.GetShortestPath(mol, acyl_carbon, ester_matches[0][-1])
        if not any(atom_idx in glycerol_atoms for atom_idx in paths):
            return False, "Ester group not connected to glycerol backbone"
    else:
        return False, "Acyl group not identified"

    # Check for remaining two hydroxyl groups or alkyl groups on the glycerol backbone
    # Count the number of hydroxyl groups attached to the glycerol carbons
    hydroxyl_count = 0
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    hydroxyl_count += 1

    if hydroxyl_count < 1 or hydroxyl_count > 2:
        return False, f"Expected 1 or 2 hydroxyl groups on glycerol, found {hydroxyl_count}"

    # Check for alkyl groups attached to glycerol carbons (excluding the esterified carbon)
    alkyl_groups = 0
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            if any(neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in glycerol_atoms for neighbor in atom.GetNeighbors()):
                alkyl_groups +=1

    total_substituents = hydroxyl_count + alkyl_groups
    if total_substituents != 2:
        return False, f"Expected total of 2 hydroxyl or alkyl groups on glycerol, found {total_substituents}"

    return True, "Contains glycerol backbone with one acyl group and two hydroxyl/alkyl groups"