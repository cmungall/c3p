"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate, consisting of a glycerol backbone
    esterified with one fatty acid chain and phosphorylated on one of the hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Glycerol backbone: Three connected carbons each with an attached oxygen (could be O or OH)
    glycerol_pattern = Chem.MolFromSmarts('OCC(O)CO')

    # Ester linkage (fatty acyl chain attached via ester bond): C(=O)OC
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C;!$(C(=O))]')

    # Phosphate group attached to glycerol: [O;H][C][O][P](=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts('COP(=O)(O)O')

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for exactly one ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, expected exactly 1"

    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Get the atom indices of glycerol backbone carbons
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    glycerol_atoms = glycerol_matches[0]
    glycerol_carbons = [atom_idx for atom_idx in glycerol_atoms if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C']

    # Check that the ester is attached to one of the glycerol carbons
    ester_match = ester_matches[0]
    ester_oxygen_idx = ester_match[2]  # Index of the oxygen atom connected to glycerol
    ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)

    ester_connected_to_glycerol = False
    for neighbor in ester_oxygen_atom.GetNeighbors():
        if neighbor.GetIdx() in glycerol_carbons:
            ester_connected_to_glycerol = True
            break
    if not ester_connected_to_glycerol:
        return False, "Ester group not attached to glycerol backbone"

    # Check that the phosphate group is attached to the glycerol backbone
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    phosphate_match = phosphate_matches[0]
    phosphate_oxygen_idx = phosphate_match[0]  # Oxygen atom connected to glycerol
    phosphate_oxygen_atom = mol.GetAtomWithIdx(phosphate_oxygen_idx)

    phosphate_connected_to_glycerol = False
    for neighbor in phosphate_oxygen_atom.GetNeighbors():
        if neighbor.GetIdx() in glycerol_carbons:
            phosphate_connected_to_glycerol = True
            break
    if not phosphate_connected_to_glycerol:
        return False, "Phosphate group not attached to glycerol backbone"

    # Check that the remaining hydroxyl group is free
    hydroxyl_groups = 0
    for carbon_idx in glycerol_carbons:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        attached_oxygens = [nbr for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == 'O']
        for oxygen in attached_oxygens:
            if oxygen.GetDegree() == 1 and oxygen.GetImplicitValence() == 1:  # Free hydroxyl
                hydroxyl_groups += 1

    if hydroxyl_groups != 1:
        return False, f"Expected one free hydroxyl group on glycerol backbone, found {hydroxyl_groups}"

    return True, "Molecule is a lysophosphatidic acid"