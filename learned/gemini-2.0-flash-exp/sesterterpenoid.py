"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are C25 isoprenoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 35 : # Increased carbon range
        return False, f"Number of carbons is {carbon_count}, expected around 25 for a sesterterpenoid"

    # Check for isoprene units (C5H8)
    isoprene_pattern = Chem.MolFromSmarts("CC(C)=CC")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    isoprene_count = len(isoprene_matches)

    # Score the presence of isoprene units.
    # Ideally, a sesterterpenoid will have 5 isoprene units (25 carbons)
    if isoprene_count < 2 : # at least 2 isoprene units
         return False, f"Too few isoprene units ({isoprene_count}). Expected at least 2."
    

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600: # Broader range
        return False, f"Molecular weight is {mol_wt}, expected between 300 and 600 for a sesterterpenoid"
    
    return True, "Matches criteria for a sesterterpenoid based on number of carbons, isoprene units, and molecular weight."