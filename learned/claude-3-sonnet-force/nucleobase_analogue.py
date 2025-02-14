"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:37585 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for heteroaromatic rings (pyrimidines, purines)
    heteroaromatic_rings = mol.GetAromaticRings(includeExocyclic=True, includeEndocyclic=True, includeEndocyclicHeteroAtoms=True, includeExoEndoCyclic=True)
    heteroaromatic_rings = [ring for ring in heteroaromatic_rings if any(atom.GetAtomicNum() in [7, 8] for atom in ring)]
    if not heteroaromatic_rings:
        return False, "No heteroaromatic rings found"
    
    # Check for common functional groups (amino, oxo, hydroxyl)
    atoms = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if 7 not in atoms or 8 not in atoms:  # Must have N and O
        return False, "Missing N or O atoms"
    
    # Check for planar structure
    calc = AllChem.MMFFGetMoleculeProperties(mol)
    if calc is None or calc.PMIBondCount < 5:  # Arbitrary threshold for planarity
        return False, "Structure is not sufficiently planar"
    
    # Check for typical nucleobase size
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    if n_heavy_atoms < 5 or n_heavy_atoms > 20:
        return False, "Size deviates too much from typical nucleobases"
    
    # If all conditions are met, classify as nucleobase analogue
    return True, "Contains heteroaromatic rings and functional groups typical of nucleobases"