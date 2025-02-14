"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.

    Minerals are typically inorganic compounds consisting of metal cations
    and inorganic anions, and are naturally occurring. This function uses
    heuristics to classify minerals based on the absence of organic carbon
    structures and the presence of metal elements.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is likely a mineral, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of metal elements (atomic numbers)
    metals = {3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,37,38,39,
              40,41,42,43,44,45,46,47,48,49,55,56,57,72,73,74,75,76,77,78,
              79,80,81,82,83,84,87,88,89,104,105,106,107,108,109,110,111,112}

    # Check for presence of metal elements
    has_metal = any(atom.GetAtomicNum() in metals for atom in mol.GetAtoms())
    if not has_metal:
        return False, "No metal elements detected"

    # Check for absence of organic carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # If carbon is bonded to hydrogen or carbon, it's likely organic
            neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            if 1 in neighbors or 6 in neighbors:
                return False, "Contains organic carbon structures"

    # Check for absence of complex organic functional groups
    organic_functional_groups = ['[CX3]=[OX1]', '[#6]#[#7]', '[#6]=[#7]']
    for fg in organic_functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_pattern):
            return False, "Contains organic functional groups"

    # Check for presence of common inorganic anions
    inorganic_anions = ['[O-]S(=O)(=O)[O-]',  # Sulfate
                        'P(=O)([O-])([O-])[O-]',  # Phosphate
                        '[C-]#N',  # Cyanide (to exclude)
                        'C(=O)([O-])[O-]',  # Carbonate
                        '[O-][N+](=O)[O-]',  # Nitrate
                        '[F-]', '[Cl-]', '[Br-]', '[I-]',  # Halides
                        '[O-]',  # Oxide
                        '[S-2]',  # Sulfide
                        '[OH-]',  # Hydroxide
                        '[O-]C(=O)[O-]'  # Oxalate (to exclude)
                        ]
    has_inorganic_anion = False
    for anion in inorganic_anions:
        anion_pattern = Chem.MolFromSmarts(anion)
        if mol.HasSubstructMatch(anion_pattern):
            has_inorganic_anion = True
            break
    if not has_inorganic_anion:
        return False, "No common inorganic anions detected"

    return True, "Likely a mineral (inorganic compound with metal elements and inorganic anions)"