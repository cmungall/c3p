"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:35741 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a benzene ring with four chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count chlorine atoms and benzene rings
    n_chlorine = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    n_benzene = len(mol.GetAromaticRings())
    
    # Tetrachlorobenzene must have exactly 4 chlorine atoms and 1 benzene ring
    if n_chlorine != 4 or n_benzene != 1:
        return False, f"Found {n_chlorine} chlorine atoms and {n_benzene} benzene rings, need exactly 4 and 1 respectively"
    
    # Check if all chlorine atoms are attached to the benzene ring
    benzene_ring = mol.GetAromaticRings()[0]
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    
    for cl_atom in chlorine_atoms:
        if not any(bond.GetBeginAtomIdx() == cl_atom.GetIdx() and bond.GetEndAtomIdx() in benzene_ring for bond in mol.GetBonds()):
            return False, "One or more chlorine atoms not attached to the benzene ring"
    
    return True, "Contains a benzene ring with 4 chlorine substituents"