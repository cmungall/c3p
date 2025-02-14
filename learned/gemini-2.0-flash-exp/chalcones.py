"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is defined as a 1,3-diphenylpropenone (benzylideneacetophenone) and its
    derivatives formed by substitution, Ar-CH=CH-C(=O)-Ar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core chalcone SMARTS pattern
    chalcone_core_pattern = Chem.MolFromSmarts("[c]~[C]=[C]~[C](=[O])~[c]")
    core_matches = mol.GetSubstructMatches(chalcone_core_pattern)

    if not core_matches:
        return False, "Core chalcone structure (Ar-CH=CH-C(=O)-Ar) not found"


    # Check molecular weight - chalcones usually > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a chalcone"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar)"