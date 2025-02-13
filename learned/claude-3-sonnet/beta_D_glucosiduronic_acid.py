"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: CHEBI:36349 beta-D-glucosiduronic acid
A glucosiduronic acid resulting from the formal condensation of any substance with beta-D-glucuronic acid to form a glycosidic bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glucuronic acid substructure
    glucuronic_acid = Chem.MolFromSmiles("O=C(O)C(O)C(O)C(O)C(O)CO")
    if not mol.HasSubstructMatch(glucuronic_acid):
        return False, "No glucuronic acid substructure found"

    # Check for glycosidic bond (acetal or ether bond to anomeric carbon)
    anomeric_carbon = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            if len(neighbors) == 2:
                anomeric_carbon = neighbors[0] if neighbors[1] == atom.GetIdx() else neighbors[1]
                break

    if anomeric_carbon is None:
        return False, "No anomeric carbon found"

    # Check that the anomeric carbon is part of a glycosidic bond
    anomeric_atom = mol.GetAtomWithIdx(anomeric_carbon)
    if anomeric_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False, "Anomeric carbon is not sp3 hybridized"

    if len(anomeric_atom.GetNeighbors()) != 2:
        return False, "Anomeric carbon does not have 2 neighbors"

    neighbor1, neighbor2 = anomeric_atom.GetNeighbors()
    if neighbor1.GetAtomicNum() != 8 or neighbor2.GetAtomicNum() != 8:
        return False, "Anomeric carbon is not part of a glycosidic bond"

    # Check that the glucuronic acid substructure is in the beta configuration
    if not AllChem.EmbedMolecule(mol):
        return False, "Could not embed molecule for 3D structure generation"

    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() == anomeric_carbon or bond.GetEndAtomIdx() == anomeric_carbon:
            begin_coords = conf.GetAtomPosition(bond.GetBeginAtomIdx())
            end_coords = conf.GetAtomPosition(bond.GetEndAtomIdx())
            v1 = begin_coords - end_coords
            v2 = end_coords - conf.GetAtomPosition(anomeric_carbon)
            cross = v1.CrossProduct(v2)
            if cross.z < 0:
                return True, "Contains beta-D-glucuronic acid substructure"

    return False, "Glucuronic acid substructure is not in the beta configuration"