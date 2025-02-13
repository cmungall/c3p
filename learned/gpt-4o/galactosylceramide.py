"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    Galactosylceramides are cerebrosides with a galactose head group and a sphingolipid backbone.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Galactose head group pattern: a hexose ring with multiple hydroxyls
    galactose_pattern = Chem.MolFromSmarts("OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"
    
    # Amide linkage pattern: -C(=O)N-
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Long-chain sphingolipid backbone pattern: -C=C-C-C-C- (with potential long aliphatic chains)
    # It's more relaxed since the sphingolipid can vary widely with chain length and saturation
    sphingolipid_pattern = Chem.MolFromSmarts("[CHX4][CHX3]([CH2])[CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(sphingolipid_pattern):
        return False, "No sphingolipid backbone found"
  
    # Optional: Count oxygen and carbon atoms to ensure they align with cerebroside structures
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if o_count < 6:
        return False, "Not enough oxygen atoms for galactose and amide functionality"
    if c_count < 30:
        return False, "Not enough carbon atoms for long sphingolipid and sugar chains"

    return True, "Contains galactose head group, amide linkage and sphingolipid backbone"