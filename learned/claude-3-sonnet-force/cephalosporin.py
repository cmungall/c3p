"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:35469 cephalosporin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are a class of beta-lactam antibiotics with a 6-membered side ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring pattern (N-C-C-N-C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[N,O]=[C-]1[N,O]C(=O)C1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Look for fused 6-membered ring (cephalosporin pattern)
    ceph_pattern = Chem.MolFromSmarts("C1[N,O]C2C(=O)N2C(=O)C1C")
    if not mol.HasSubstructMatch(ceph_pattern):
        return False, "No 6-membered cephalosporin ring found"

    # Check for other common cephalosporin features
    thiazole_pattern = Chem.MolFromSmarts("c1nccs1")
    tetrazole_pattern = Chem.MolFromSmarts("c1cn[nH]n1")
    
    has_thiazole = mol.HasSubstructMatch(thiazole_pattern)
    has_tetrazole = mol.HasSubstructMatch(tetrazole_pattern)
    
    if has_thiazole or has_tetrazole:
        return True, "Contains beta-lactam and 6-membered cephalosporin ring, with thiazole or tetrazole group"

    # Additional checks based on molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, "Molecular weight outside typical cephalosporin range"
    
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Insufficient rotatable bonds for cephalosporin"

    return True, "Contains beta-lactam and 6-membered cephalosporin ring with appropriate molecular properties"