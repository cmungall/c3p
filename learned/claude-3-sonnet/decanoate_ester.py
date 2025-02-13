"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:33673 decanoate ester

A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for decanoate group (decanoyl, CCCCCCCCCC(=O)-) attached to an oxygen
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)[O;X2]")
    decanoate_matches = mol.GetSubstructMatches(decanoate_pattern)
    if not decanoate_matches:
        return False, "No decanoate group found"

    # Check if the oxygen is part of an alcohol or phenol group
    alcohol_pattern = Chem.MolFromSmarts("[O;X2][C;X4]")
    phenol_pattern = Chem.MolFromSmarts("[O;X2][c]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    if not alcohol_matches and not phenol_matches:
        return False, "Decanoate group not attached to an alcohol or phenol"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for decanoate esters"

    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2 or oxygen_count > 8:
        return False, "Number of oxygens outside typical range for decanoate esters"

    return True, "Contains a decanoate group attached to an alcohol or phenol group"