"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:33839 thiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a -SH group attached to an aliphatic or aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfhydryl (-SH) group
    sh_pattern = Chem.MolFromSmarts("[SH]")
    sh_matches = mol.GetSubstructMatches(sh_pattern)
    if not sh_matches:
        return False, "No sulfhydryl (-SH) group found"
    
    # Check if sulfhydryl is attached to aliphatic or aromatic carbon
    is_thiol = False
    for atom_idx in sh_matches[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetTotalNumHs() == 1:  # Check for -SH
            neighbors = atom.GetNeighbors()
            if any(nb.GetAtomicNum() == 6 and nb.GetHybridization() in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP3) for nb in neighbors):
                is_thiol = True
                break
    
    if is_thiol:
        return True, "Contains sulfhydryl (-SH) group attached to aliphatic or aromatic carbon"
    else:
        return False, "Sulfhydryl (-SH) group not attached to aliphatic or aromatic carbon"