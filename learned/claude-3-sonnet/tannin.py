"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:18035 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_tannin(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds containing galloyl or hexahydroxydiphenoyl (HHDP) moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for galloyl or HHDP substructures
    galloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1C(=O)O")
    hhdp_pattern = Chem.MolFromSmarts("O=C1OC2=C(O)C(O)=C(O)C(O)=C2C(O)=C1O")
    if not mol.HasSubstructMatch(galloyl_pattern) and not mol.HasSubstructMatch(hhdp_pattern):
        return False, "No galloyl or HHDP moieties found"
    
    # Check for polyphenolic nature
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if aromatic_rings < 2 or hydroxyl_groups < 5:
        return False, "Not sufficiently polyphenolic"
    
    # Check molecular weight range (typically 500-3000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 3000:
        return False, "Molecular weight outside typical tannin range"
    
    # Additional checks based on specific tannin classes
    # (e.g., procyanidins, proanthocyanidins, ellagitannins)
    # ...

    return True, "Contains galloyl or HHDP moieties and exhibits polyphenolic nature"