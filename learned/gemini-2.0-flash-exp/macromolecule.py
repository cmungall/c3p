"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:25367 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a large molecule built from repeating smaller units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Check number of heavy (non-H) atoms
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    # Calculate number of rotatable bonds - indicates chain length
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Calculate number of bonds (faster proxy for chain length)
    num_bonds = mol.GetNumBonds()

    # calculate number of rings to discard cyclic molecules
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Apply heuristics for macromolecule classification
    if mol_wt > 1800:
        return True, "Very high molecular weight, likely a macromolecule."
    
    if num_heavy_atoms > 120 or num_rotatable_bonds > 40 or num_bonds > 150:
          return True, "High molecular weight or large number of heavy atoms or rotatable bonds or number of bonds, indicating a macromolecule."
    
    elif num_heavy_atoms > 150 and num_rotatable_bonds > 45:
          return True, "Large number of heavy atoms and rotatable bonds, indicating a macromolecule"
    elif num_rotatable_bonds > 70 or num_bonds > 200 :
           return True, "High number of rotatable bonds or bonds, indicating a macromolecule."

    elif num_rings > 10:
          return False, "Too many rings, not likely a macromolecule."
    
    else:
        return False, "Does not meet macromolecule criteria based on molecular weight, heavy atoms, number of rotatable bonds, number of bonds, or rings."