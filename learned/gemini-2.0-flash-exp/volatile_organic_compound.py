"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound based on its SMILES string.
    This is an approximation based on molecular properties and may not be 100% accurate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate properties
    mol_wt = Descriptors.MolWt(mol)
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    num_hbd = Descriptors.NumHDonors(mol)
    num_hba = Descriptors.NumHAcceptors(mol)

    # Heuristic-based rules
    if mol_wt > 300:
         return False, "Molecular weight too high, likely not a VOC"
    
    if num_heavy_atoms > 20:
        return False, "Too many heavy atoms, likely not a VOC"
    
    if (num_hbd + num_hba) > 3:
         return False, "Too many hydrogen bond donors and acceptors, likely not a VOC"
    
    return True, "Meets criteria for a likely VOC"