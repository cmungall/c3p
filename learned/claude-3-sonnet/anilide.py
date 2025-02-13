"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:35759 anilide
Anilides are aromatic amides obtained by acylation of aniline.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    # Look for aromatic ring attached to amide nitrogen
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    for amide_match in amide_matches:
        n_atom = mol.GetAtomWithIdx(amide_match[1])
        if any(mol.GetAtomWithIdx(neighbor).IsAromatic() for neighbor in n_atom.GetNeighbors()):
            break
    else:
        return False, "No aromatic ring attached to amide nitrogen"

    # Check for additional acyl group attached to amide nitrogen
    acyl_pattern = Chem.MolFromSmarts("C(=O)")
    for amide_match in amide_matches:
        n_atom = mol.GetAtomWithIdx(amide_match[1])
        if any(mol.GetAtomWithIdx(neighbor).HasSubstructMatch(acyl_pattern)
               for neighbor in n_atom.GetNeighbors()):
            break
    else:
        return False, "No acyl group attached to amide nitrogen"

    return True, "Contains an aromatic amide group derived from acylation of aniline"