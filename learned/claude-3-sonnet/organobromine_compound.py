"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:33636 Organobromine compounds
An organobromine compound is a compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of bromine atoms
    bromine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35]
    if not bromine_atoms:
        return False, "No bromine atoms found"
    
    # Check for carbon-bromine bonds
    has_c_br_bond = any(atom.GetTotalNumHs() == 0 and
                        any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors())
                        for atom in bromine_atoms)
    
    if has_c_br_bond:
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bonds found"