"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are characterized by a 30-carbon skeleton, often with multiple fused rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a triterpenoid, False otherwise
         str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a minimum number of rings using rdMolDescriptors
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
      return False, f"Too few rings, got {num_rings}, need at least 4."

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 35:
       return False, f"Carbon count out of range {c_count}. Triterpenoids must have ~30 C."
        
    # Check oxygen count
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 0 or o_count > 15:
        return False, f"Oxygen count out of range {o_count}, requires between 0 to 15 oxygens."


    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 1000:
        return False, f"Molecular weight out of range {mol_wt}, must be between 350 and 1000."


    return True, "Meets criteria for a triterpenoid (multiple rings, appropriate C count, O count, and molecular weight)"