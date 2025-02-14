"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: CHEBI:36973 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
import numpy as np

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
    This function checks for the presence of a repeated monosaccharide pattern.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get adjacency matrix
    adj_matrix = GetAdjacencyMatrix(mol)
    
    # Check for repeating patterns
    n_atoms = mol.GetNumAtoms()
    patterns = []
    for start in range(n_atoms):
        for length in range(3, n_atoms//2 + 1):  # Minimum length 3 to avoid small fragments
            pattern = [int(x) for x in list(adj_matrix[start, start:start+length])]
            if pattern not in patterns:
                patterns.append(pattern)
            else:
                # Found repeating pattern, check if it is a monosaccharide
                submol = Chem.PathToSubmol(mol, np.array(pattern), atomMap=list(range(start, start+length)))
                if submol.GetNumHeavyAtoms() > 3 and len(submol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H][CX4][CX3](O)[CX3](O)[CX3](O)[CX3][CX3]"))) > 0:
                    # Repeating pattern matches monosaccharide
                    return True, f"Found repeating monosaccharide pattern of length {length}"
    
    return False, "No repeating monosaccharide pattern found"