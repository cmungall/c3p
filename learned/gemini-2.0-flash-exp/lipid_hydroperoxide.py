"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid containing one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for hydroperoxy group (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX2H1]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy group found"

    # 2. Check for a long carbon chain (at least 8 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
         return False, "Too few carbons for a lipid"
    
    # 3. Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
         return False, "Molecular weight too low for a lipid"

    return True, "Contains at least one hydroperoxy group and a long carbon chain"