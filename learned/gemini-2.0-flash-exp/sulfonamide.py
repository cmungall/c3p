"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide has the structure R-S(=O)2-N-R'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for the sulfonamide core
    sulfonamide_core_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    sulfonamide_nitrogen_pattern = Chem.MolFromSmarts("S(=O)(=O)[N;!O]")
    
    # Search for substructure match
    if not mol.HasSubstructMatch(sulfonamide_core_pattern):
         return False, "Sulfonamide core not found"

    if not mol.HasSubstructMatch(sulfonamide_nitrogen_pattern):
      return False, "Sulfonamide nitrogen is bonded to oxygen only"
      
    #Verify nitrogen is attached to at least one C or H
    matches = mol.GetSubstructMatches(sulfonamide_core_pattern)
    for match in matches:
        nitrogen_atom_idx = match[1]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_atom_idx)
        n_neighbor_is_carbon_or_hydrogen = False
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() in [1,6]:
              n_neighbor_is_carbon_or_hydrogen = True
              break
        if not n_neighbor_is_carbon_or_hydrogen:
            amide_pattern = Chem.MolFromSmarts("[N;!H0](C=O)")
            if not mol.HasSubstructMatch(amide_pattern):
              return False, "Nitrogen not bonded to at least one C or H and not an amide"


    return True, "Contains a sulfonamide group (S(=O)2-N) with valid nitrogen"