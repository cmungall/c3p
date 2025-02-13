"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:36567 tetraterpenoid

A tetraterpenoid is any terpenoid derived from a tetraterpene. The term includes compounds in which
the C40 skeleton of the parent tetraterpene has been rearranged or modified by the removal of one
or more skeletal atoms (generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate number of atoms
    num_atoms = mol.GetNumAtoms()
    
    # Tetraterpenoids typically have 40-50 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 40 or num_carbons > 50:
        return False, f"Number of carbon atoms ({num_carbons}) outside expected range for tetraterpenoid"
    
    # Check for long carbon chains and rings, typical of terpenoids
    Sssr = Chem.GetSymmSSSR(mol)
    has_long_chains = any(len(ring) > 8 for ring in Sssr)
    has_rings = len(Sssr) > 0
    if not has_long_chains or not has_rings:
        return False, "No long carbon chains or rings found, typical of terpenoids"
    
    # Check for isoprene units (C5H8)
    isoprene_pattern = Chem.MolFromSmarts("[C@H]([CH3])=C[CH2][CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 8:  # At least 8 isoprene units expected for tetraterpenoids
        return False, f"Found {len(isoprene_matches)} isoprene units, expected at least 8"
    
    # Check molecular weight - tetraterpenoids typically 500-700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 700:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside expected range for tetraterpenoid"
    
    return True, "Contains typical features of a tetraterpenoid"