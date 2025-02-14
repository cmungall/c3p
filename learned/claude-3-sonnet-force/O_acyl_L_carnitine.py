"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: CHEBI:36584 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is a carnitine molecule with an acyl group attached to the
    oxygen atom, and the carnitine component must have the L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carnitine backbone pattern (C-C-C-N-C)
    carnitine_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][N+]([C])[C]")
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "No carnitine backbone found"

    # Look for acyl group attached to oxygen (-O-C(=O)-)
    acyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, expected exactly 1"

    # Check for L-configuration of carnitine
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnspecified=True)
    if not chiral_centers:
        return False, "No chiral center found for carnitine"

    chiral_atom = mol.GetAtomWithIdx(chiral_centers[0][0])
    if chiral_atom.GetProp("_CIPCode") != "S":
        return False, "Carnitine is not in L-configuration"

    return True, "Contains carnitine backbone with an acyl group attached to oxygen, L-configuration"