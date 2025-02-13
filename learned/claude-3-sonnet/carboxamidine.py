"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35813 carboxamidine

Carboxamidines are compounds having the structure RC(=NR)NR2. The term is used 
as a suffix in systematic nomenclature to denote the -C(=NH)NH2 group including 
its carbon atom.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxamidine pattern: R-C(=NR)-NR2
    carboxamidine_pattern = Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]")
    if not mol.HasSubstructMatch(carboxamidine_pattern):
        return False, "No carboxamidine substructure found"
    
    # Check if all atoms are accounted for in the matched pattern
    matched_atoms = set()
    for match in mol.GetSubstructMatches(carboxamidine_pattern):
        matched_atoms.update(match)
    if len(matched_atoms) != mol.GetNumAtoms():
        return False, "Carboxamidine substructure does not cover the entire molecule"
    
    return True, "Molecule contains the carboxamidine substructure R-C(=NR)-NR2"