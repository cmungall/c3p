"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the glycerophospholipid backbone
    # Glycerol backbone with two fatty acid chains esterified at sn-1 and sn-2 positions, and phosphate at sn-3
    glycerophospholipid_pattern = Chem.MolFromSmarts("""
    [$([C@@H](CO[P](=O)(O)O)(O[C](=O)[C;H1,$(C([CH2])[CH2])][C;H1,$(C([CH2])[CH2])[CH2])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]),
    $([C@@H](CO[P](=O)(O)O)(O[C](=O)[C;H1,$(C([CH2])[CH2])][C;H1,$(C([CH2])[CH2])[CH2])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2])])]
    """)

    if not mol.HasSubstructMatch(glycerophospholipid_pattern):
        return False, "No glycerophospholipid backbone with fatty acids and phosphate group found"

    # Define SMARTS pattern for the inositol ring attached via phosphate group
    # Inositol ring: six-membered ring with carbons each bearing a hydroxyl group
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check if inositol ring is connected via phosphate group to glycerol backbone
    # We can check for phosphatidylinositol substructure
    phosphatidylinositol_pattern = Chem.MolFromSmarts("[O]-[P](=O)([O])-[O]-C1(CO)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(phosphatidylinositol_pattern):
        return False, "No phosphatidylinositol structure found"

    # Find the inositol ring atoms
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"
    inositol_atoms = inositol_matches[0]
    inositol_atom_set = set(inositol_atoms)

    # Count the number of phosphate groups attached to the inositol ring hydroxyls
    num_phosphates_on_inositol = 0

    for atom_idx in inositol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look for hydroxyl groups on inositol ring carbons
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if oxygen is bonded to phosphorus
                for nei in neighbor.GetNeighbors():
                    if nei.GetAtomicNum() == 15 and nei.GetIdx() not in inositol_atom_set:
                        num_phosphates_on_inositol += 1
                        break

    # The linkage phosphate to glycerol counts as one, so additional phosphates indicate phosphorylation
    if num_phosphates_on_inositol >= 2:
        num_additional_phosphates = num_phosphates_on_inositol - 1
        return True, f"Phosphoinositide with {num_additional_phosphates} additional phosphate group(s) on inositol ring"
    else:
        return False, "No additional phosphate groups on inositol ring"