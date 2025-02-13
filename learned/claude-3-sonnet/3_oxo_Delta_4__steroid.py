"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:51106 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid with a 3-oxo group and a C=C double bond
    at the alpha,beta position (between C4 and C5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("C(=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "No 3-oxo group found"
    elif len(oxo_matches) > 1:
        return False, "More than one oxo group found"
    
    # Check for alpha,beta C=C double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    alpha_beta_match = False
    for match in double_bond_matches:
        bond = mol.GetBondBetweenAtoms(match[0], match[1])
        if bond.GetBeginAtomIdx() == 3 and bond.GetEndAtomIdx() == 4:
            alpha_beta_match = True
            break
    if not alpha_beta_match:
        return False, "No alpha,beta C=C double bond found"
    
    # Additional checks
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1:
        return False, "Molecule is too rigid for a steroid"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, "Molecular weight outside typical range for steroids"
    
    return True, "Contains steroid backbone with 3-oxo group and alpha,beta C=C double bond"