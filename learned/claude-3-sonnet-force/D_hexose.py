"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16646 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

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

    # Check if molecule contains exactly 6 carbon atoms and 6 oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count != 6:
        return False, "Not a hexose (does not contain 6 carbons and 6 oxygens)"

    # Check for the presence of a ring structure
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "Not a hexose (no ring structure found)"

    # Get the stereochemistry of the molecule
    stereo_chem = Chem.FindMolChiralCenters(mol)

    # Define a reference D-hexose (e.g., D-glucose) and its stereochemistry
    reference_smiles = "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"
    reference_mol = Chem.MolFromSmiles(reference_smiles)
    reference_stereo_chem = Chem.FindMolChiralCenters(reference_mol)

    # Compare the stereochemistry of the molecule with the reference D-hexose
    if stereo_chem == reference_stereo_chem:
        return True, "Molecule has the same stereochemistry as a reference D-hexose"
    else:
        return False, "Stereochemistry does not match a reference D-hexose"

# Example usage
print(is_D_hexose("OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"))  # True, 'Molecule has the same stereochemistry as a reference D-hexose'
print(is_D_hexose("OC[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O"))  # False, 'Stereochemistry does not match a reference D-hexose'